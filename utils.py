def read_data(path) :
    points = []

    f = open(path, "r")
    data_string = f.read()

    for point in data_string.split('{') :
        if len(point) > 0 :
            point_split = point.split(',')
            point = [float(point_split[0]), float(point_split[1]), float(point_split[2].split('}')[0])]

            points.append(point)

    return points



if __name__ == '__main__':
    path = './data/sensors.txt'

    data = read_data(path)
    print(len(data))
    print(data)