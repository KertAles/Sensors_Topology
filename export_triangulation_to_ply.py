def export_ply(triangulation, points):
    triangulation = [t for t in triangulation
                     if len(t) == 3]
    data = """ply
format ascii 1.0           
comment Exported triangulation
comment 
element vertex {}
property float x
property float y
property float z
element face {}
property list uchar int vertex_index
end_header
""".format(len(points), len(triangulation))
    for point in points:
        data += " ".join(map(str, point)) + "\n"
    for simplex in triangulation:
        line = "{} ".format(len(simplex))
        line += " ".join(map(str, simplex))
        line += "\n"
        data += line
    return data    
