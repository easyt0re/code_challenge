import math

def cross_prod(coord_a, coord_b):
	# given 2 vectors on the same plane
	# the cross product is the plane norm vector
    result = [coord_a[1]*coord_b[2] - coord_a[2]*coord_b[1],
            coord_a[2]*coord_b[0] - coord_a[0]*coord_b[2],
            coord_a[0]*coord_b[1] - coord_a[1]*coord_b[0]]

    return result

def dot_prod(coord_a, coord_b):
	result = coord_a[0] * coord_b[0] + coord_a[1] * coord_b[1] + coord_a[2] * coord_b[2]
	# result = math.sumprod(coord_a, coord_b)

	return result

def add_list(coord_a, coord_b):
    result =    [
                    coord_a[0] + coord_b[0],
                    coord_a[1] + coord_b[1],
                    coord_a[2] + coord_b[2]
                ]
    
    return result

def scale_list(coord_a, scalar_K):
    result =    [
                    coord_a[0] * scalar_K,
                    coord_a[1] * scalar_K,
                    coord_a[2] * scalar_K
                ]

    return result

def diff_list(coord_a, coord_b):
    result = add_list(coord_a, scale_list(coord_b, -1.0))

    return result

def unit_norm_vector(coord_a, coord_b):
	norm_vector = cross_prod(coord_a, coord_b)
	result = norm_vector / (dot_prod(norm_vector, norm_vector)) ** 0.5

	return result

def norm_list(coord_a):
    return (dot_prod(coord_a, coord_a)) ** 0.5

# math.dist
# def pt2plane(coord_a, coord_b, u_norm_vector):
# 	# use this to compute height
# 	dist_height = dot_prod(coord_b - coord_a, u_norm_vector)
# 	proj_pt = coord_b - dist_height * u_norm_vector

# 	return (proj_pt, dist_height)

class dist_btw_2_coords(object):
    """docstring for dist_btw_2_coords"""
    def __init__(self, coord_a, coord_b):
        super(dist_btw_2_coords, self).__init__()
        self.coord_a = coord_a
        self.coord_b = coord_b
        self.dist_btw = math.dist(coord_a, coord_b)

# class abs_vect(object):
#     """docstring for abs_vect"""
#     def __init__(self, coord_a, **kwds):
#         super().__init__()
#         self.coord_a = coord_a
#         print("am i here")

#     def add_vect(self, vect_a):
#         pass

#     def scaling(self, scalar_K):
#         pass
        
# class rel_vect(abs_vect):
#     """docstring for rel_vect"""
#     def __init__(self, coord_b, **kwds):
#         super().__init__()
#         print('im here')
#         # self.coord_a = coord_a
#         self.coord_b = coord_b

#         self.dist_btw = math.dist(coord_a, coord_b)

#     def dot_prod(self):
#         pass

#     def cross_prod(self):
#         pass

class simple_vect(object):
    """docstring for simple_vect"""
    def __init__(self, coord_a, coord_b=[0, 0, 0]):
        super().__init__()
        self.coord_a = coord_a
        self.coord_b = coord_b

        self.dist_btw = math.dist(coord_a, coord_b)

    def add_vect(self, vect_a):
        pass

    def scaling(self, scalar_K):
        pass

    def dot_prod(self):
        pass

    def cross_prod(self):
        pass

def find_circle(coord_list_3Pts):
    [coord_a, coord_b, coord_c] = coord_list_3Pts
    vect_ab = simple_vect(coord_a, coord_b)
    vect_ac = simple_vect(coord_a, coord_c)
    vect_cb = simple_vect(coord_c, coord_b)

    if (vect_ab.dist_btw > vect_ac.dist_btw) and (vect_ab.dist_btw > vect_cb.dist_btw):
        largest_dist_vect = vect_ab
        coord_third = coord_list_3Pts[2]
    elif vect_ac.dist_btw > vect_cb.dist_btw: #if reached this point, bc or ac must be >= ab
        largest_dist_vect = vect_ac
        coord_third = coord_list_3Pts[1]
    else:
        largest_dist_vect = vect_cb 
        coord_third = coord_list_3Pts[0]

    # largest_dist_vect = find_max_dist(vect_ab, vect_ac, vect_cb)
    coord_mid = scale_list(add_list(largest_dist_vect.coord_a, largest_dist_vect.coord_b), 1/2)
    if math.dist(coord_mid, coord_third) > largest_dist_vect.dist_btw / 2:
        coord_center, circle_r = fine_DIY(coord_list_3Pts)
    else:
        coord_center = coord_mid
        circle_r = largest_dist_vect.dist_btw / 2

    return coord_center, circle_r, unit_norm_vector

# def find_max_dist(vect_ab, vect_ac, vect_cb):
#     if (vect_ab.dist_btw > vect_ac.dist_btw) and (vect_ab.dist_btw > vect_cb.dist_btw):
#         return vect_ab
#     elif vect_ac.dist_btw > vect_cb.dist_btw: #if reached this point, bc or ac must be >= ab
#         return vect_ac
#     else:
#         return vect_cb

def fine_DIY(coord_list_3Pts):
    # math here: https://math.stackexchange.com/a/1743505
    p1, p2, p3 = coord_list_3Pts
    diff_21 = diff_list(p2, p1)
    diff_31 = diff_list(p3, p1)
    norm_vect = cross_prod(diff_21, diff_31)
    u_x_axis = scale_list(diff_21, 1 / norm_list(diff_21))
    u_z_axis = scale_list(norm_vect, 1 / norm_list(norm_vect))
    u_y_axis = cross_prod(u_z_axis, u_x_axis)
    # TODO: some repeated code here

    diff_21_px = dot_prod(diff_21, u_x_axis)
    diff_31_px = dot_prod(diff_31, u_x_axis)
    diff_31_py = dot_prod(diff_31, u_y_axis)

    center_px = diff_21_px / 2
    center_py = ((diff_31_px - center_px) ** 2 + diff_31_py ** 2 - center_px ** 2) / 2 / diff_31_py

    coord_center = add_list(p1, add_list(scale_list(u_x_axis, center_px), scale_list(u_y_axis, center_py)))
    radius_2d = ((diff_31_px - center_px) ** 2 + (diff_31_py - center_py) ** 2) ** 0.5
    radius_3d = math.dist(coord_center, p1)

    return coord_center, [radius_2d, radius_3d]

class simple_plane(object):
    """docstring for simple_plane"""
    def __init__(self, coord_list_3Pts):
        super().__init__()
        p1, p2, p3 = coord_list_3Pts
        diff_21 = diff_list(p2, p1)
        diff_31 = diff_list(p3, p1)
        diff_32 = diff_list(p3, p2)
        norm_vect = cross_prod(diff_21, diff_31)
        u_x_axis = scale_list(diff_21, 1 / norm_list(diff_21))
        u_z_axis = scale_list(norm_vect, 1 / norm_list(norm_vect))
        u_y_axis = cross_prod(u_z_axis, u_x_axis)
        self.origin = p1

        self.u_x_axis = u_x_axis
        self.u_z_axis = u_z_axis
        self.u_y_axis = u_y_axis
        self.diff_21 = diff_21
        self.diff_31 = diff_31
        self.diff_32 = diff_32

        # find circle(s)
        # find 3p circle
        diff_21_px = dot_prod(diff_21, u_x_axis)
        diff_31_px = dot_prod(diff_31, u_x_axis)
        diff_31_py = dot_prod(diff_31, u_y_axis)

        center_px = diff_21_px / 2
        center_py = ((diff_31_px - center_px) ** 2 + diff_31_py ** 2 - center_px ** 2) / 2 / diff_31_py

        coord_center = add_list(p1, add_list(scale_list(u_x_axis, center_px), scale_list(u_y_axis, center_py)))
        # radius_2d = ((diff_31_px - center_px) ** 2 + (diff_31_py - center_py) ** 2) ** 0.5
        radius_3d = math.dist(coord_center, p1)
        self.circle3p = simple_circle(coord_center, radius_3d)

        # find 2p circle
        vect_ab = simple_vect(p1, p2)
        vect_ac = simple_vect(p1, p3)
        vect_cb = simple_vect(p3, p2)

        if (vect_ab.dist_btw > vect_ac.dist_btw) and (vect_ab.dist_btw > vect_cb.dist_btw):
            largest_dist_vect = vect_ab
            coord_third = coord_list_3Pts[2]
        elif vect_ac.dist_btw > vect_cb.dist_btw: #if reached this point, bc or ac must be >= ab
            largest_dist_vect = vect_ac
            coord_third = coord_list_3Pts[1]
        else:
            largest_dist_vect = vect_cb 
            coord_third = coord_list_3Pts[0]

        coord_mid = scale_list(add_list(largest_dist_vect.coord_a, largest_dist_vect.coord_b), 1/2)
        if math.dist(coord_mid, coord_third) > largest_dist_vect.dist_btw / 2:
            self.hasCircle2p = False
        else:
            self.hasCircle2p = True
            self.circle2p = simple_circle(coord_mid, largest_dist_vect.dist_btw / 2)

class simple_circle(object):
    """docstring for simple_circle"""
    def __init__(self, coord_center, circle_r):
        super(simple_circle, self).__init__()
        self.coord_center = coord_center
        self.circle_r = circle_r
        


def proj_3d_to_2d(planeIn, starPtOut):
    vector_to_plane = [starPtOut[0] - planeIn.origin[0], starPtOut[1] - planeIn.origin[1], starPtOut[2] - planeIn.origin[2]]
    x_on_plane = dot_prod(vector_to_plane, planeIn.u_x_axis)
    y_on_plane = dot_prod(vector_to_plane, planeIn.u_y_axis)
    height = dot_prod(planeIn.u_z_axis, vector_to_plane)
    # positive, on the norm vector side
    # starPtOnPlane = vector_to_plane - height * planeIn.u_norm_vector

    return ((x_on_plane, y_on_plane), height)

def height_and_dist(coord_center, circle_r, coord_another, unit_norm_vector):
    vect2plane = diff_list(coord_another, coord_center)
    height = dot_prod(vect2plane, unit_norm_vector)
    dist_vect = diff_list(vect2plane, scale_list(unit_norm_vector, height))
    distance = norm_list(dist_vect)

    return height, distance

def check_circle(a_circle, unit_norm_vector, coord_another):
    coord_center = a_circle.coord_center
    circle_r = a_circle.circle_r
    
    vect2plane = diff_list(coord_another, coord_center)
    height = dot_prod(vect2plane, unit_norm_vector)
    dist_vect = diff_list(vect2plane, scale_list(unit_norm_vector, height))
    distance = norm_list(dist_vect)

    return height, distance



# Test with three points in 3D space
p1 = [1, 0, 0]
p2 = [0, 2, 0]
p3 = [0, 0, 3]

coord_1 = [0, 0, 0]
coord_2 = [0, 0, 1]
coord_3 = [1, 1, 0]
coord_4 = [1, 1, 1]
coord_5 = [0, 1, 0]
# print("input done")
# print(add_list(p2, p3))
# print(scale_list(p2, -3))
# print(dot_prod(p2, p3))
# print(diff_list(p2, p1))
# print(diff_list(p3, p1))
# print(norm_list(p3))
# print(simple_vect(p1).coord_a)
# print(simple_vect(p1, p2).dist_btw)

# center, radius = fine_DIY(p1, p2, p3)
# print(f"Center: {center}, Radius: {radius}")

coord_list = [coord_1, coord_2, coord_3, coord_4, coord_5]

# test_case_1
coord_1 = [1, 0, 0]
coord_2 = [1, 1, 0]
coord_3 = [0, 0, 0]
coord_4 = [0, 0, 1]
test_case_1 = [coord_1, coord_2, coord_3, coord_4]
# result = 1.57079633

# test_case_2
coord_1 = [-100, 0, 0]
coord_2 = [10, 0, 10]
coord_3 = [-10, -10, -10]
coord_4 = [0, 0, 0]
test_case_2 = [coord_1, coord_2, coord_3, coord_4]
# result = 41938.65135885

# test_case_1
coord_1 = [10, 20, 30]
coord_2 = [0, 0, 0]
coord_3 = [-100, 1000, -20]
coord_4 = [100, -20, 33]
coord_5 = [8, -7, 900]
coord_6 = [-100, -223, -23]
coord_7 = [3, 0, 3]
test_case_3 = [coord_1, coord_2, coord_3, coord_4, coord_5, coord_6, coord_7]
# result = 298192571.11934924

test_case_in = test_case_3

min_volumn = 1e100
from itertools import combinations
# import bisect

for each_comb in combinations(test_case_in, 3):
    # this is for every combination for 3 out of N points

    # plane as class
    planeIn = simple_plane(each_comb)
    # # rely on func
    # coord_center, circle_r, unit_norm_vector = find_circle(each_comb)

    # always check 3p circle
    for each_coord in coord_list:
        # plane as class
        height, distance = check_circle(planeIn.circle3p, planeIn.u_z_axis, each_coord)
        # # rely on func
        # height, distance = height_and_dist(coord_center, circle_r, each_coord, unit_norm_vector)
        if height != 0:
            if 'prev_height' in locals():
                # this is not the first point
                if height * prev_height < 0:
                    # this is not a bottom circle, break out
                    break
                else:
                    if distance > planeIn.circle3p.circle_r:
                        # this is not in circle, break out
                        break
                    else:
                        # put height into a list
                        # maybe try to keep a sorted list
                        # bisect.insort(height_list, height)

                        # actually, only the largest is enough
                        new_height = abs(height)
                        if new_height > max_height:
                            max_height = new_height
                        
            else:
                # this is the first point
                max_height = abs(height)

            prev_height = height

        else:
            # this is either the 3 points or 4 points co-planar
            pass

    # end of for loop 3p height search 
    if 'max_height' in locals():
        volumn = math.pi * planeIn.circle3p.circle_r ** 2 * max_height
        if volumn < min_volumn:
            min_volumn = volumn

        # end of 3p circle check, reset things
        # intentionally put in if structure to make sure they exist
        del prev_height, max_height

    # conditionally check 2p circle, same as 3p
    if planeIn.hasCircle2p:
        for each_coord in coord_list:
            # plane as class
            height, distance = check_circle(planeIn.circle2p, planeIn.u_z_axis, each_coord)
            # # rely on func
            # height, distance = height_and_dist(coord_center, circle_r, each_coord, unit_norm_vector)
            if height != 0:
                if 'prev_height' in locals():
                    # this is not the first point
                    if height * prev_height < 0:
                        # this is not a bottom circle, break out
                        break
                    else:
                        if distance > planeIn.circle2p.circle_r:
                            # this is not in circle, break out
                            break
                        else:
                            # put height into a list
                            # maybe try to keep a sorted list
                            # bisect.insort(height_list, height)

                            # actually, only the largest is enough
                            new_height = abs(height)
                            if new_height > max_height:
                                max_height = new_height
                            
                else:
                    # this is the first point
                    max_height = abs(height)

                # update "memory"
                prev_height = height

            else:
                # this is either the 3 points or 4 points co-planar
                pass

        # end of for loop 2p height search 
        if 'max_height' in locals():
            volumn = math.pi * planeIn.circle2p.circle_r ** 2 * max_height
            if volumn < min_volumn:
                min_volumn = volumn

            # end of 2p circle check, reset things
            # intentionally put in if structure to make sure they exist
            del prev_height, max_height

print(min_volumn)