context("Test declutter function")


test_that("After decluttering (+sorting) by any number of points, if any decluttered coordinates are still one apart", {

    ## A few examples
    coords.to.check = list(c(8,10,15,16,19), c(20,8,9,19,18))

    ## see if any are one apart /after/ decluttering.
    for(coords in coords.to.check){
        processed.coords = declutter(coords, how.close = 1, sort=T, indexonly = F)
        if(length(processed.coords) >= 2){
            one.apart.exists = any(abs(processed.coords[1:(length(processed.coords)-1)]
                                       - processed.coords[2:length(processed.coords)])==1)
        }
        expect_equal(one.apart.exists, FALSE)
    }
}
