#include <Rcpp.h>

#include <iostream>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef std::vector< std::vector<unsigned> > FamilyRowMap;

static double sumBranches(unsigned ndx, double product, const Rcpp::NumericVector &sharingProbs,
const double observedPValue);

static bool contains(const Rcpp::CharacterVector &vec, const std::string &str)
{
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        if (Rcpp::as<std::string>(vec[i]) == str)
        {
            return true;
        }
    }
    return false;
}

static void checkInputs(const Rcpp::IntegerMatrix &alleles,
const Rcpp::CharacterVector &variants, const Rcpp::CharacterVector &famIds,
const Rcpp::NumericVector &sharingProbs, const std::string &filter,
const Rcpp::IntegerVector &minorAllele, double alpha)
{
    if (variants.size() != alleles.ncol())
    {
        Rcpp::stop("Mismatch between variant names and the SnpMatrix");
    }
    if (famIds.size() != alleles.nrow())
    {
        Rcpp::stop("Mistmatch between family ids and the SnpMatrix");
    }
    if (!filter.empty() && (alpha <= 0.0 || alpha > 1.0))
    {
        Rcpp::stop("Invalid alpha value for chosen filter");
    }
    for (unsigned i = 0; i < famIds.size(); ++i)
    {
        if (!contains(sharingProbs.attr("names"), Rcpp::as<std::string>(famIds[i])))
        {
            Rcpp::stop("Can't find family ID in the sharing probabilities");
        }
    }
}

// This maps each family id to all of the rows in the SnpMatrix that contain an affected
// subject of that family id. The map is indexed based on the sharingProbs vector,
// i.e. if you want to find the rows corresponding to the subjects in sharingProbs[i] you
// can use familyRowMap[i]
static FamilyRowMap getFamilyRowMap(const Rcpp::NumericVector &sharingProbs,
const Rcpp::CharacterVector &famIds)
{
    FamilyRowMap familyRowMap;
    Rcpp::CharacterVector uniqueFamIds = sharingProbs.attr("names");
    for (unsigned i = 0; i < uniqueFamIds.size(); ++i)
    {
        familyRowMap.push_back(std::vector<unsigned>());
        std::string thisFamId = Rcpp::as<std::string>(uniqueFamIds[i]);
        for (unsigned j = 0; j < famIds.size(); ++j) // loop over rows of SnpMatrix
        {
            if (famIds[j] == thisFamId)
            {
                familyRowMap[i].push_back(j);
            }
        }
    }
    return familyRowMap;
}

static bool minorAllelePresent(const Rcpp::IntegerMatrix &snpMat, int minorAllele,
const std::vector<unsigned> &rowNdx, unsigned colNdx)
{
    for (unsigned i = 0; i < rowNdx.size(); ++i)
    {
        if (snpMat(rowNdx[i], colNdx) == 2 || snpMat(rowNdx[i], colNdx) == minorAllele)
        {
            return true;
        }
    }
    return false;
}

static bool minorAlleleShared(const Rcpp::IntegerMatrix &snpMat, int minorAllele,
const std::vector<unsigned> &rowNdx, unsigned colNdx)
{
    for (unsigned i = 0; i < rowNdx.size(); ++i)
    {
        if (snpMat(rowNdx[i], colNdx) != minorAllele && snpMat(rowNdx[i], colNdx) != 2)
        {
            return false;
        }
    }
    return true;
}

static Rcpp::NumericVector calculatePotentialPValues(const Rcpp::IntegerMatrix &snpMat,
const Rcpp::CharacterVector &famIds, const Rcpp::NumericVector &sharingProbs,
const Rcpp::IntegerVector &minorAllele, const FamilyRowMap &familyRowMap,
const std::vector<unsigned> &validVariants)
{
    // calculate potential p-values
    unsigned nVariants = static_cast<unsigned>(validVariants.size());
    Rcpp::NumericVector ppvals(nVariants, 1.0);
    #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
    for (unsigned i = 0; i < nVariants; ++i)
    {
        unsigned varNdx = validVariants[i];
        for (unsigned j = 0; j < sharingProbs.size(); ++j)
        {
            if (minorAllelePresent(snpMat, minorAllele[varNdx], familyRowMap[j], varNdx))
            {
                ppvals[i] *= sharingProbs[j];
            }
        }
    }
    return ppvals;
}

static double findPValueCutoff(const Rcpp::NumericVector &ppvals, const std::string &filter, double alpha)
{
    if (filter.empty())
    {
        return 1.0;
    }
    Rcpp::NumericVector sorted = Rcpp::clone(ppvals);
    sorted.sort();
    unsigned ndx = 0;
    while (sorted[ndx] < alpha / static_cast<double>(ndx + 1))
    {
        ++ndx;
    }
    return sorted[ndx - 1];
}

static Rcpp::NumericVector filterVector(const Rcpp::NumericVector &vecIn, double minVal)
{
    Rcpp::NumericVector vecOut;
    for (unsigned i = 0; i < vecIn.size(); ++i)
    {
        if (vecIn[i] >= minVal)
        {
            vecOut.push_back(vecIn[i]);
        }
    }
    return vecOut;
}

struct ProductIndexPair
{
    double product;
    unsigned index;
    ProductIndexPair(double p, unsigned n) : product(p), index(n) {}
};

static double sumBranches_nr(unsigned ndx, double product, const Rcpp::NumericVector &sharingProbs,
const double observedPValue)
{
    double pvalue = 0.0;
    std::vector<ProductIndexPair> stack;
    std::vector<double> signifProbs;
    stack.push_back(ProductIndexPair(product, ndx));
    while (!stack.empty())
    {
        ProductIndexPair current = stack.back();
        stack.pop_back();
        if (current.index < sharingProbs.size())
        {
            double leftProduct = sharingProbs[current.index] * current.product;
            double rightProduct = (1.0 - sharingProbs[current.index]) * current.product;
            if (leftProduct <= observedPValue)
            {
                signifProbs.push_back(leftProduct);
                pvalue += leftProduct;
            }
            else
            {
                stack.push_back(ProductIndexPair(leftProduct, current.index + 1));
            }
            if (rightProduct <= observedPValue)
            {
                signifProbs.push_back(rightProduct);
                pvalue += rightProduct;
            }
            else
            {
                stack.push_back(ProductIndexPair(rightProduct, current.index + 1));
            }
        }
    }
    double pval2 = 0.0;
    for (unsigned i = 0; i < signifProbs.size(); ++i)
    {
        //printf("%.12f, ", signifProbs[i]);
        pval2 += signifProbs[i];
    }
    //printf("\n");
    //if (sharingProbs.size() > 1)
    //{
    //    unsigned whichMaxSharingProb = 0;
    //    double maxSharingProb = sharingProbs[0];
    //    for (unsigned i = 0; i < sharingProbs.size(); ++i)
    //    {
    //        if (sharingProbs[i] > maxSharingProb)
    //        {
    //            maxSharingProb = sharingProbs[i];
    //            whichMaxSharingProb = i;
    //        }
    //    }
    //    double oneFalse = 1.0;
    //    for (unsigned i = 0; i < sharingProbs.size(); ++i)
    //    {
    //        if (i != whichMaxSharingProb)
    //        {
    //            oneFalse *= sharingProbs[i];
    //        }
    //        else
    //        {
    //            //oneFalse *= (1.0 - sharingProbs[i]);
    //            oneFalse *= sharingProbs[i];
    //        }
    //    }
    //    //std::cout << sharingProbs.size() << " Observed: " << observedPValue
    //    //    << ", One False: " << oneFalse << '\n';
    //}
    //if (signifProbs.size() > 1)
    //{
    //    //std::cout << "nValues: " << signifProbs.size() << ", pval1: " << pvalue << ", pval2: " << pval2 << '\n';
    //}
    return pvalue;
}

static double sumBranches(unsigned ndx, double product, const Rcpp::NumericVector &sharingProbs,
const double observedPValue)
{
    if (ndx >= sharingProbs.size())
    {
        return 0;
    }
    double leftProduct = sharingProbs[ndx] * product;
    double rightProduct = (1.0 - sharingProbs[ndx]) * product;
    double leftSum = (leftProduct <= observedPValue + 0.001)
        ? leftProduct
        : sumBranches(ndx + 1, leftProduct, sharingProbs, observedPValue);
    double rightSum = (rightProduct <= observedPValue + 0.001)
        ? rightProduct
        : sumBranches(ndx + 1, rightProduct, sharingProbs, observedPValue);
    return leftSum + rightSum;
}

// [[Rcpp::export]]
double multipleFamilyPValue_cpp(const Rcpp::NumericVector &sharingProbs,
const Rcpp::LogicalVector &observedSharing, double minPValue=0.0)
{
    double observedPValue = 1.0;
    for (unsigned i = 0; i < observedSharing.size(); ++i)
    {
        observedPValue *= observedSharing[i] ? sharingProbs[i] : (1.0 - sharingProbs[i]);
    }
    if (observedPValue == 0.0)
    {
        return 0.0;
    }
    double pvalue = 1.0;
    double burnInObservedPValue = (minPValue > 0.0) ? 1.0 : observedPValue * 2.0;
    while (pvalue > minPValue && burnInObservedPValue > observedPValue)
    {
        burnInObservedPValue = std::max(burnInObservedPValue / 10.0, observedPValue);
        pvalue = sumBranches_nr(0, 1, sharingProbs, burnInObservedPValue);
    }
    return pvalue;
}

static Rcpp::NumericVector calculatePValues(const Rcpp::IntegerMatrix &snpMat,
const Rcpp::CharacterVector &famIds, 
const Rcpp::NumericVector &sharingProbs, const Rcpp::IntegerVector &minorAllele,
double cutoff, const FamilyRowMap &familyRowMap, const Rcpp::NumericVector &ppvals,
const std::vector<unsigned> &validVariants)
{
    unsigned nVariants = static_cast<unsigned>(validVariants.size());
    Rcpp::NumericVector pvals(nVariants, -1.0);
    for (unsigned i = 0; i < nVariants; ++i)
    {
        unsigned varNdx = validVariants[i];
        if (ppvals[i] <= cutoff)
        {
            Rcpp::NumericVector probs;
            Rcpp::LogicalVector pattern;
            for (unsigned j = 0; j < familyRowMap.size(); ++j)
            {
                if (minorAllelePresent(snpMat, minorAllele[varNdx], familyRowMap[j], varNdx))
                {
                    probs.push_back(sharingProbs[j]);
                    pattern.push_back(minorAlleleShared(snpMat, minorAllele[varNdx], familyRowMap[j], varNdx));
                }
            }
            pvals[i] = multipleFamilyPValue_cpp(probs, pattern);
        }
    }
    return pvals;
}

static Rcpp::IntegerVector determineMinorAlleles(const Rcpp::IntegerMatrix &snpMat)
{
    Rcpp::IntegerVector minorAllele;
    for (int j = 0; j < snpMat.ncol(); ++j)
    {
        int count = 0;
        for (int i = 0; i < snpMat.nrow(); ++i)
        {
            if (snpMat(i,j) == 1)
                ++count;
            else if (snpMat(i,j) == 3)
                --count;
        }
        int minorAlleleCandidate = (count > 0) ? 3 : 1;
        minorAllele.push_back(0);
        for (int i = 0; i < snpMat.nrow(); ++i)
        {
            if (snpMat(i,j) == minorAlleleCandidate || snpMat(i,j) == 2)
            {
                minorAllele[j] = minorAlleleCandidate;
            }
        }
    }
    return minorAllele;
}

static std::vector<unsigned> findValidVariants(const Rcpp::IntegerVector &minorAllele)
{
    std::vector<unsigned> validNdx;
    for (unsigned j = 0; j < minorAllele.size(); ++j)
    {
        if (minorAllele[j] != 0)
        {
            validNdx.push_back(j);
        }
    }
    return validNdx;
}

static Rcpp::CharacterVector getValidVariantnames(const std::vector<unsigned> &validVariants,
const Rcpp::CharacterVector &allVariantNames)
{
    Rcpp::CharacterVector validNames;
    for (unsigned i = 0; i < validVariants.size(); ++i)
    {
        validNames.push_back(allVariantNames[validVariants[i]]);
    }
    return validNames;
}

// [[Rcpp::export]]
Rcpp::List multipleVariantPValue_cpp(const Rcpp::IntegerMatrix &snpMat,
const Rcpp::CharacterVector &variants, const Rcpp::CharacterVector &famIds,
const Rcpp::NumericVector &sharingProbs,
const Rcpp::Nullable<Rcpp::CharacterVector> &rfilter,
const Rcpp::Nullable<Rcpp::IntegerVector> &minorAlleleInput,
double alpha)
{
#ifdef _OPENMP
    std::cout << "Running on " << omp_get_max_threads() << " threads\n";
#endif
    // determine minor allele value for each variant
    Rcpp::IntegerVector minorAllele = minorAlleleInput.isNull()
        ? determineMinorAlleles(snpMat)
        : Rcpp::IntegerVector(minorAlleleInput);
    std::vector<unsigned> validVariants = findValidVariants(minorAllele);
    unsigned nIgnoredVariants = snpMat.ncol() - validVariants.size();
    if (nIgnoredVariants > 0)
    {
        std::cout << "Ignoring " << nIgnoredVariants << " variants not present in any subject\n";
    }

    // parse arguments from R and check validity of parameters
    std::string filter = (rfilter.isNull()) ? "" : Rcpp::as<std::string>(rfilter);
    checkInputs(snpMat, variants, famIds, sharingProbs, filter, minorAllele, alpha);

    // calculate p-values for each variant
    FamilyRowMap familyRowMap(getFamilyRowMap(sharingProbs, famIds));
    Rcpp::NumericVector potPValues = calculatePotentialPValues(snpMat, famIds,
        sharingProbs, minorAllele, familyRowMap, validVariants);
    double cutoff = findPValueCutoff(potPValues, filter, alpha);
    Rcpp::NumericVector pvalues = calculatePValues(snpMat, famIds, sharingProbs,
        minorAllele, cutoff, familyRowMap, potPValues, validVariants);
    
    // return p-values along with indices of variants used
    return Rcpp::List::create(
        Rcpp::Named("pvalues") = pvalues,
        Rcpp::Named("potential_pvalues") = potPValues,
        Rcpp::Named("valid_variants") = Rcpp::wrap(validVariants)
    );
}

// [[Rcpp::export]]
double enrichmentPValue_cpp(const Rcpp::IntegerMatrix &snpMat,
const Rcpp::CharacterVector &famIds, const Rcpp::NumericVector &sharingProbs,
const Rcpp::IntegerVector &minorAllele, double threshold)
{
    Rcpp::NumericVector probs;
    Rcpp::LogicalVector pattern;
    unsigned nVariants = static_cast<unsigned>(snpMat.ncol());
    FamilyRowMap familyRowMap(getFamilyRowMap(sharingProbs, famIds));
    for (unsigned i = 0; i < nVariants; ++i)
    {
        for (unsigned j = 0; j < sharingProbs.size(); ++j)
        {
            if (minorAllelePresent(snpMat, minorAllele[i], familyRowMap[j], i))
            {
                probs.push_back(sharingProbs[j]);
                pattern.push_back(minorAlleleShared(snpMat, minorAllele[i], familyRowMap[j], i));
            }
        }
    }
    return multipleFamilyPValue_cpp(probs, pattern, threshold);
}
