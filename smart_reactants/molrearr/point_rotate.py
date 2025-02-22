import numpy as np

# 1、绕Z周旋转gamma角

def rotate_Z(x, y, z,  gamma):
    gamma = gamma * (np.pi / 180)
    x_r = np.cos(gamma)*x - np.sin(gamma)*y
    y_r = np.sin(gamma)*x + np.cos(gamma)*y
    z_r = z
    print(f"{(x, y, z)} rotate {gamma*(180/np.pi)} degrees around the Z-axis,result {(x_r, y_r, z_r)}")
    return x_r, y_r, z_r


# 2、绕Y轴旋转beta角
def rotate_Y(x, y, z, beta):
    beta = beta * (np.pi / 180)
    x_r = np.cos(beta)*x + np.sin(beta)*z
    y_r = y
    z_r = -np.sin(beta)*x + np.cos(beta)*z
    print(f"{(x, y, z)} rotate {beta*(180/np.pi)} degrees around the Y-axis,result {(x_r, y_r, z_r)}")
    return x_r, y_r, z_r


# 3、绕X轴旋转alpha角
def rotate_X(x, y, z, alpha):
    alpha = alpha * (np.pi / 180)
    x_r = x
    y_r = np.cos(alpha)*y - np.sin(alpha)*z
    z_r = np.sin(alpha)*y + np.cos(alpha)*z
    print(f"{(x, y, z)} rotate {alpha*(180/np.pi)} degrees around the X-axis,result {(x_r, y_r, z_r)}")
    return x_r, y_r, z_r



# https://zh.wikipedia.org/wiki/%E6%97%8B%E8%BD%AC%E7%9F%A9%E9%98%B5
# 可能是坐标系的不同，维基百科上的斜对角的sin的负号，和我们使用的刚好是交换过来的

if __name__ == '__main__':
    # 1、绕Z轴旋转，Z值不变
    # 点(2, 0, 12)绕Z轴旋转90度，旋转后 (0, 2, 12)
    rotate_Z(2, 0, 12, 90)

    # 2、绕Y轴旋转，Y值不变
    # 点(3, 4, 0)绕Y轴顺时针旋转90度，旋转后 (0, 4, -3)
    rotate_Y(3, 4, 0, 90)

    # 3、绕X轴旋转，X值不变
    # 点(5, 0, 7)绕X轴顺时针旋转90度，旋转后 (5, -7, 0)
    rotate_X(5, 0, 7, 90)

    x, y, z = 1, 2, 3
    x, y, z = rotate_X(x, y, z, 90)
    x, y, z = rotate_Y(x, y, z, -90)
    x, y, z = rotate_Z(x, y, z, -90)
    # 结果：(1,2,3) -> (-3, 2, 1)

    x, y, z = 1, 2, 3
    x, y, z = rotate_Z(x, y, z, 90)
    x, y, z = rotate_Y(x, y, z, -90)
    x, y, z = rotate_X(x, y, z, -90)
    # 结果： (1,2,3) - > (-3, -2, -1)
	# 因此不同的旋转顺序会导致旋转的结果不一样
