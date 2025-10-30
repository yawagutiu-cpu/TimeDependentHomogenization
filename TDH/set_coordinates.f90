subroutine set_coordinates

    use set_parameter
    
    implicit none
    
    gauss(1, 1:3) = (/-ga, -ga, -ga/)
    gauss(2, 1:3) = (/ ga, -ga, -ga/)
    gauss(3, 1:3) = (/ ga,  ga, -ga/)
    gauss(4, 1:3) = (/-ga,  ga, -ga/)
    gauss(5, 1:3) = (/-ga, -ga,  ga/)
    gauss(6, 1:3) = (/ ga, -ga,  ga/)
    gauss(7, 1:3) = (/ ga,  ga,  ga/)
    gauss(8, 1:3) = (/-ga,  ga,  ga/)

    polar(1, 1:3) = (/-1.0d0, -1.0d0, -1.0d0/)
    polar(2, 1:3) = (/ 1.0d0, -1.0d0, -1.0d0/)
    polar(3, 1:3) = (/ 1.0d0,  1.0d0, -1.0d0/)
    polar(4, 1:3) = (/-1.0d0,  1.0d0, -1.0d0/)
    polar(5, 1:3) = (/-1.0d0, -1.0d0,  1.0d0/)
    polar(6, 1:3) = (/ 1.0d0, -1.0d0,  1.0d0/)
    polar(7, 1:3) = (/ 1.0d0,  1.0d0,  1.0d0/)
    polar(8, 1:3) = (/-1.0d0,  1.0d0,  1.0d0/)

end subroutine