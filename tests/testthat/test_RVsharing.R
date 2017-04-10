context('Top Level Testing')

test_that('test old RVsharing implementation',
{
    data(testPed)
    expect_equal(RVsharing(testPedList[[1]])@pshare, 1/31)
    expect_equal(RVsharing(testPedList[[2]])@pshare, 1/63)
    expect_equal(RVsharing(testPedList[[3]])@pshare, 1/255)
    expect_equal(RVsharing(testPedList[[4]])@pshare, 1/127)
    expect_equal(RVsharing(testPedList[[5]])@pshare, 0.00004096178)
    expect_equal(RVsharing(testPedList[[6]])@pshare, 1/745)
    expect_equal(RVsharing(testPedList[[7]])@pshare, 1/31)
    expect_equal(unname(RVsharing(testPedList[[8]])@pshare), 1/7)
    expect_equal(RVsharing(testPedList[[11]])@pshare, 1/255)
    expect_equal(RVsharing(testPedList[[12]])@pshare, 1/361)
    expect_equal(RVsharing(testPedList[[13]])@pshare, 1/127)
    expect_equal(RVsharing(testPedList[[15]])@pshare, 1/31)
    expect_equal(RVsharing(testPedList[[16]])@pshare, 1/31)

    expect_error(RVsharing(testPedList[[9]]))
})
