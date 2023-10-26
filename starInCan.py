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
        super(simple_vect, self).__init__()
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
    # coord_list_3Pts = [coord_a, coord_b, coord_c]

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
    print("axis")
    print(u_x_axis)
    print(u_y_axis)
    print(u_z_axis)

    diff_21_px = dot_prod(diff_21, u_x_axis)
    diff_31_px = dot_prod(diff_31, u_x_axis)
    diff_31_py = dot_prod(diff_31, u_y_axis)

    center_px = diff_21_px / 2
    center_py = ((diff_31_px - center_px) ** 2 + diff_31_py ** 2 - center_px ** 2) / 2 / diff_31_py

    coord_center = add_list(p1, add_list(scale_list(u_x_axis, center_px), scale_list(u_y_axis, center_py)))
    radius_2d = ((diff_31_px - center_px) ** 2 + (diff_31_py - center_py) ** 2) ** 0.5
    radius_3d = math.dist(coord_center, p1)

    return coord_center, [radius_2d, radius_3d]

class Plane:
 def __init__(self, starPt1, starPt2, starPt3):

   self.origin = starPt1
   vector_a = [starPt2[0] - starPt1[0], starPt2[1] - starPt1[1], starPt2[2] - starPt1[2]]
   vector_b = [starPt3[0] - starPt1[0], starPt3[1] - starPt1[1], starPt3[2] - starPt1[2]]
   self.u_x_axis = get_unit_vector(vector_a)
   self.u_z_axis = get_unit_vector(cross_prod(vector_a, vector_b))
   self.u_y_axis = get_unit_vector(cross_prod(self.u_z_axis, self.u_x_axis))

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



# Test with three points in 3D space
p1 = [1, 0, 0]
p2 = [0, 2, 0]
p3 = [0, 0, 3]

coord_1 = [0, 0, 0]
coord_2 = [0, 0, 1]
coord_3 = [1, 1, 0]
coord_4 = [1, 1, 1]
coord_5 = [0, 1, 0]
print("input done")
print(add_list(p2, p3))
print(scale_list(p2, -3))
print(dot_prod(p2, p3))
print(diff_list(p2, p1))
print(diff_list(p3, p1))
print(norm_list(p3))
print(simple_vect(p1).coord_a)
print(simple_vect(p1, p2).dist_btw)

center, radius = fine_DIY(p1, p2, p3)
print(f"Center: {center}, Radius: {radius}")

coord_list = [coord_1, coord_2, coord_3, coord_4, coord_5]
min_volumn = 1e10
from itertools import combinations
for each_comb in combinations(coord_list, 3):
    coord_center, circle_r, unit_norm_vector = find_circle(each_comb)
    for each_coord in coord_list:
        height, distance = height_and_dist(coord_center, circle_r, each_coord, unit_norm_vector)
        # TODO: pseudo code from this point on
        if height != 0:
            if height * prev_height < 0:
                # this is not a bottom circle, break out
                pass
            else:
                if distance > circle_r:
                    # this is not in circle, break out
                    pass
                else:
                    # put height into a list
                    # maybe try to keep a sorted list
                    pass

    # end of height search
    if height_list is not empty:
        volumn = pi * circle_r ** 2 * max(height_list)
        if volumn < min_volumn:
            min_volumn = volumn


# # by chatGPT

# def calculate_circle_np(p1, p2, p3):
#     # Convert points to numpy arrays
#     p1, p2, p3 = np.array(p1), np.array(p2), np.array(p3)

#     # Calculate the cross product of vectors
#     cp = np.cross(p1-p2, p1-p3)

#     # Calculate the square of magnitudes
#     d = 2 * np.sum(cp**2)

#     # Calculate the center of the circle
#     center = (np.sum((p1**2)*(p2-p3) + (p2**2)*(p3-p1) + (p3**2)*(p1-p2)) / d) * cp

#     # Calculate the radius of the circle
#     radius = np.sqrt(np.sum((p1 - center)**2))

#     return center.tolist(), radius

# # # Test with three points in 3D space
# # p1 = [0, 0, 0]
# # p2 = [1, 0, 0]
# # p3 = [0, 1, 0]

# # center, radius = calculate_circle_center_radius(p1, p2, p3)
# # print(f"Center: {center}, Radius: {radius}")


# def calculate_circle_center_radius(p1, p2, p3):
#     ax, ay, az = p1
#     bx, by, bz = p2
#     cx, cy, cz = p3

#     dax = (ax - cx)
#     day = (ay - cy)
#     daz = (az - cz)

#     dbx = (bx - cx)
#     dby = (by - cy)
#     dbz = (bz - cz)

#     d = 2.0 * (dax * (dby * daz - dbz * day) - day * (dbx * daz - dbz * dax) + daz * (dbx * day - dby * dax))
#     diff_12 = diff_list(p1, p2)
#     diff_23 = diff_list(p2, p3)
#     diff_31 = diff_list(p3, p1)
#     norm_vect = cross_prod(diff_12, scale_list(diff_31, -1.0))
#     vect_mag = dot_prod(norm_vect, norm_vect)
#     d = 2 * vect_mag
#     print(d)
#     # part1 = scale_list(diff_23, dot_prod(p1, p1))
#     # part2 = scale_list(diff_31, dot_prod(p2, p2))
#     # part3 = scale_list(diff_12, dot_prod(p3, p3))
#     # parts_together = add_list(add_list(part1, part2), part3)
#     # coord_center = scale_list(norm_vect, math.fsum(parts_together) / d)

#     ux = ((ax * ax + ay * ay + az * az) * (dby * daz - dbz * day) + 
#           (bx * bx + by * by + bz * bz) * (day * dbz - daz * dby) +
#           (cx * cx + cy * cy + cz * cz) * (daz * dby - day * dbz)) / d

#     uy = ((ax * ax + ay * ay + az * az) * (dbz * dax - dbx * daz) +
#           (bx * bx + by * by + bz * bz) * (daz * dax - day * dbz) +
#           (cx * cx + cy * cy + cz * cz) * (day * dbx - dax * dby)) / d

#     uz = ((ax * ax + ay * ay + az * az) * (dbx * day - dby* dax) +
#           (bx* bx + by* by + bz* bz) *(dax* day - dbx* daz) +
#           (cx* cx + cy* cy + cz* cz) *(dby* daz - day* dbz)) / d

#     radius = ((ux-ax)**2+(uy-ay)**2+(uz-az)**2)**0.5
#     return [ux, uy, uz], radius
#     # radius = math.dist(p1, coord_center)
#     # return coord_center, radius