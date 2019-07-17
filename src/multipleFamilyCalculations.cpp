#include <Rcpp.h>

#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef std::vector< std::vector<unsigned> > FamilyRowMap;

static double leafSum(unsigned ndx, double product, const Rcpp::NumericVector &sharingProbs,
const double observedPValue);

static double sumBranches(unsigned ndx, double product, const Rcpp::NumericVector &sharingProbs,
const double observedPValue);

bool contains(const Rcpp::CharacterVector &vec, const std::string &str)
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

void checkInputs(const Rcpp::IntegerMatrix &alleles,
const Rcpp::CharacterVector &variants, const Rcpp::CharacterVector &famIds,
const Rcpp::NumericVector &sharingProbs, const std::string &filter,
const Rcpp::NumericVector &minorAllele, double alpha)
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
FamilyRowMap getFamilyRowMap(const Rcpp::NumericVector &sharingProbs,
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

Rcpp::NumericVector calculatePotentialPValues(const Rcpp::IntegerMatrix &snpMat,
const Rcpp::CharacterVector &famIds, const Rcpp::NumericVector &sharingProbs,
const Rcpp::NumericVector &minorAllele, const FamilyRowMap &familyRowMap)
{
    // calculate potential p-values
    unsigned nVariants = static_cast<unsigned>(snpMat.ncol());
    Rcpp::NumericVector ppvals(nVariants, 1.0);
    #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static)
    for (unsigned i = 0; i < nVariants; ++i)
    {
        for (unsigned j = 0; j < sharingProbs.size(); ++j)
        {
            if (minorAllelePresent(snpMat, minorAllele[i], familyRowMap[j], i))
            {
                ppvals[i] *= sharingProbs[j];
            }
        }
    }
    return ppvals;
}

double findPValueCutoff(const Rcpp::NumericVector &ppvals, const std::string &filter, double alpha)
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

static double leafSum(unsigned ndx, double product, const Rcpp::NumericVector &sharingProbs,
const double observedPValue)
{
    return (product < observedPValue + 0.001)
        ? product
        : sumBranches(ndx + 1, product, sharingProbs, observedPValue);
}

static double sumBranches(unsigned ndx, double product, const Rcpp::NumericVector &sharingProbs,
const double observedPValue)
{
    if (ndx >= sharingProbs.size())
    {
        return 0;
    }
    double leftProduct = sharingProbs[ndx] * product;
    double rightProduct = (1 - sharingProbs[ndx]) * product;
    return leafSum(ndx, leftProduct, sharingProbs, observedPValue)
        + leafSum(ndx, rightProduct, sharingProbs, observedPValue);
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
    double burnInObservedPValue = 1.0;
    while (pvalue > minPValue && burnInObservedPValue > observedPValue)
    {
        burnInObservedPValue = std::max(burnInObservedPValue / 10.0, observedPValue);
        pvalue = sumBranches(0, 1, sharingProbs, burnInObservedPValue);
    }
    return pvalue;
}

Rcpp::NumericVector calculatePValues(const Rcpp::IntegerMatrix &snpMat,
const Rcpp::CharacterVector &famIds, 
const Rcpp::NumericVector &sharingProbs, const Rcpp::NumericVector &minorAllele,
double cutoff, const FamilyRowMap &familyRowMap, const Rcpp::NumericVector &ppvals)
{
    unsigned nVariants = static_cast<unsigned>(snpMat.ncol());
    Rcpp::NumericVector pvals(nVariants, -1.0);
    for (unsigned i = 0; i < nVariants; ++i)
    {
        if (ppvals[i] <= cutoff)
        {
            Rcpp::NumericVector probs;
            Rcpp::LogicalVector pattern;
            for (unsigned j = 0; j < familyRowMap.size(); ++j)
            {
                if (minorAllelePresent(snpMat, minorAllele[i], familyRowMap[j], i))
                {
                    probs.push_back(sharingProbs[j]);
                    pattern.push_back(minorAlleleShared(snpMat, minorAllele[i], familyRowMap[j], i));
                }
            }
            pvals[i] = multipleFamilyPValue_cpp(probs, pattern);
        }
    }
    return pvals;
}

// [[Rcpp::export]]
Rcpp::List multipleVariantPValue_cpp(const Rcpp::IntegerMatrix &alleles,
const Rcpp::CharacterVector &variants, const Rcpp::CharacterVector &famIds,
const Rcpp::NumericVector &sharingProbs,
const Rcpp::Nullable<Rcpp::CharacterVector> &rfilter,
const Rcpp::NumericVector &minorAllele,
double alpha)
{
    std::string filter = (rfilter.isNull()) ? "" : Rcpp::as<std::string>(rfilter);
    checkInputs(alleles, variants, famIds, sharingProbs, filter, minorAllele, alpha);
    FamilyRowMap familyRowMap(getFamilyRowMap(sharingProbs, famIds));
    Rcpp::NumericVector potPValues = calculatePotentialPValues(alleles, famIds,
        sharingProbs, minorAllele, familyRowMap);
    double cutoff = findPValueCutoff(potPValues, filter, alpha);
    Rcpp::NumericVector pvalues = calculatePValues(alleles, famIds, sharingProbs,
        minorAllele, cutoff, familyRowMap, potPValues);
    pvalues = filterVector(pvalues, 0.0);
    return Rcpp::List::create(
        Rcpp::Named("pvalues") = pvalues,
        Rcpp::Named("potential_pvalues") = potPValues
    );
}

// [[Rcpp::export]]
double enrichmentPValue_cpp(const Rcpp::IntegerMatrix &snpMat,
const Rcpp::CharacterVector &famIds, const Rcpp::NumericVector &sharingProbs,
const Rcpp::NumericVector &minorAllele, double threshold)
{
    Rcpp::NumericVector probs;
    Rcpp::LogicalVector pattern;
    FamilyRowMap familyRowMap(getFamilyRowMap(sharingProbs, famIds));
    for (unsigned i = 0; i < snpMat.ncol(); ++i)
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
