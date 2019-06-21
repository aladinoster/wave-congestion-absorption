import os
import numpy as np
from lxml import etree as ET


def RotMat(theta):
    return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])


def RotateVector(p, theta):
    return RotMat(theta) @ p


def RotationsVector(p, n_max=10, k=4):
    l_rot = []
    for n in range(0, n_max):
        theta = 2 * np.pi * n / (n_max * k)
        l_rot.append(RotateVector(p, theta))
    return l_rot


def create_internal_points(points: list, center: tuple, n_points: int = 10, k=4) -> tuple:
    """Compute a set of internal points for a ring-road 
    
    
    :return: list of points within the circle
    :rtype: tuple 
    """

    vct_pts = [np.array(pt) for pt in points]
    center = np.array(center)
    vct_ctr = [pt - center for pt in vct_pts]
    int_points = [RotationsVector(v, n_points, k)[1:] + center for v in vct_ctr]

    return int_points


def extract_element_troncons(filename):
    tree = ET.parse(filename)
    root = tree.getroot()
    branch_tree = "RESEAUX/RESEAU/TRONCONS"
    troncons = root.xpath(branch_tree)[0].getchildren()
    return troncons


def extract_points(filename):
    troncons = extract_element_troncons(filename)
    pointstxt = [tron.attrib.get("extremite_amont") for tron in troncons[1:-1]]
    points = [tuple(float(v) for v in ptxt.split(" ")) for ptxt in pointstxt]
    return points


def inject_xml_data(filename):
    troncons = extract_element_troncons(filename)[1:-1]
    points = extract_points(filename)
    center = np.add.reduce(points) / 4
    int_points = create_internal_points(points, center)
    for i, tp in enumerate(zip(int_points, troncons)):
        xml_file = create_xml_data(tp[0])
        tp[1].append(xml_file)
    return troncons


def create_xml_data(intern_points) -> ET.ElementTree:
    """[summary]
    
    :param x_intern_coor: [description]
    :type x_intern_coor: [type]
    :param y_intern_coor: [description]
    :type y_intern_coor: [type]
    :return: No return 
    :rtype: None
    """
    root = ET.Element("POINTS_INTERNES")

    for x, y in intern_points:
        ET.SubElement(root, "POINT_INTERNE", coordonnees=f"{x} {y}")

    tree = ET.ElementTree(root)
    return tree
    # tree.write("filename.xml")


def create_xml_random() -> ET.ElementTree:
    # import xml.etree.ElementTree as ET

    # create the file structure
    data = ET.Element("data")
    items = ET.SubElement(data, "items")
    item1 = ET.SubElement(items, "item")
    item2 = ET.SubElement(items, "item")
    item1.set("name", "item1")
    item2.set("name", "item2")
    item1.text = "item1abc"
    item2.text = "item2abc"

    # create a new XML file with the results
    # mydata = ET.tostring(data)
    # myfile = open("items2.xml", "w")
    # myfile.write(mydata)

    tree = ET.ElementTree(data)
    return tree


def fix_values(element, dct_values) -> ET.ElementTree:
    for k, v in dct_values.items():
        element.set(k, v)
    return element


def create_xml_simulation() -> ET.ElementTree:

    # Create a general scheme
    attr_qname = ET.QName("http://www.w3.org/2001/XMLSchema-instance", "noNamespaceSchemaLocation")
    nsmap = {"xsi": "http://www.w3.org/2001/XMLSchema-instance"}
    root = ET.Element("ROOT_SYMUBRUIT", {attr_qname: "reseau.xsd"}, version="2.05", nsmap=nsmap)

    # Fill simulaiton
    simulations = ET.SubElement(root, "SIMULATIONS")
    simulation = ET.SubElement(simulations, "SIMULATION")
    DCT_SIMULATION = {
        "id": "simID",
        "pasdetemps": "0.1",
        "debut": "07:00:00",
        "fin": "07:15:00",
        "loipoursuite": "exacte",
        "comportementflux": "iti",
        "date": "2019-06-19",
        "titre": "ring_simulation",
        "proc_deceleration": "false",
        "seed": "1",
    }
    simulation = fix_values(simulation, DCT_SIMULATION)

    # Add traffic
    trafics = ET.SubElement(root, "TRAFICS")
    trafic = ET.SubElement(trafics, "TRAFIC")
    DCT_TRAFIC = {"id": "trafID", "accbornee": "true", "coeffrelax": "4", "chgtvoie_ghost": "false"}
    trafic = fix_values(trafic, DCT_TRAFIC)

    # Add troncon
    troncons = ET.SubElement(trafic, "TRONCONS")
    TP_TRONCONS = ("Zone_A", "Zone_B", "Zone_C", "Zone_D", "Zone_E", "Zone_F")
    for tron in TP_TRONCONS:
        troncon = ET.SubElement(troncons, "TRONCON")
        troncon.set("id", tron)

    # Add vehicle_type
    vehtypes = ET.SubElement(trafic, "TYPES_DE_VEHICULE")
    TP_VEHTYPES = (
        {"id": "CAV", "w": "-5", "kx": "0.12", "vx": "25"},
        {"id": "HDV", "w": "-5", "kx": "0.12", "vx": "25"},
    )
    TP_ACCEL = (
        {"ax": "1.5", "vit_sup": "5.8"},
        {"ax": "1", "vit_sup": "8"},
        {"ax": "0.5", "vit_sup": "infini"},
    )
    for veh in TP_VEHTYPES:
        veh_type = ET.SubElement(vehtypes, "TYPE_DE_VEHICULE")
        veh_type = fix_values(veh_type, veh)
        acc_maps = ET.SubElement(veh_type, "ACCELERATION_PLAGES")
        # Add acceleration
        for acc in TP_ACCEL:
            acc_map = ET.SubElement(acc_maps, "ACCELERATION_PLAGE")
            acc_map = fix_values(acc_map, acc)

    # Adding borders
    extremes = ET.SubElement(trafic, "EXTREMITES")
    extreme = ET.SubElement(extremes, "EXTREMITE")
    TP_EXTREMES = ("Ext_In", "Ext_Out")
    for ext in TP_EXTREMES:
        extreme.set("id", ext)
        if ext == "Ext_In":
            flux_global = ET.SubElement(extreme, "FLUX_GLOBAL")
            flux = ET.SubElement(flux_global, "FLUX")

    # Adding internal connections
    connections = ET.SubElement(trafic, "CONNEXIONS_INTERNES")
    TP_CONNECTIONS = ("A_to_B", "B_to_C", "C_to_D", "D_to_E")
    for con in TP_CONNECTIONS:
        connection = ET.SubElement(connections, "CONNEXION_INTERNE")
        connection.set("id", con)

    # Adding network links + internal points
    networks = ET.SubElement(root, "RESEAUX")
    network = ET.SubElement(networks, "RESEAU")
    DCT_NETWORK_INFO = {"id": "resID"}

    network = fix_values(network, DCT_NETWORK_INFO)

    radious = 200

    TP_LINKS = (
        {
            "id": "Zone_A",
            "id_eltamont": "Ext_In",
            "id_eltaval": "A_to_B",
            "extremite_amont": f"-{np.pi * radious} 0",
            "extremite_aval": "0 0",
        },
        {
            "id": "Zone_B",
            "id_eltamont": "A_to_B",
            "id_eltaval": "B_to_C",
            "extremite_amont": "0 0",
            "extremite_aval": f"{radious} {radious}",
        },
        {
            "id": "Zone_C",
            "id_eltamont": "B_to_C",
            "id_eltaval": "C_to_D",
            "extremite_amont": f"{radious} {radious}",
            "extremite_aval": f"0 {2 * radious}",
        },
        {
            "id": "Zone_D",
            "id_eltamont": "C_to_D",
            "id_eltaval": "D_to_E",
            "extremite_amont": f"0 {2 * radious}",
            "extremite_aval": f"{-radious} {radious}",
        },
        {
            "id": "Zone_E",
            "id_eltamont": "D_to_E",
            "id_eltaval": "A_to_B",
            "extremite_amont": f"{-radious} {radious}",
            "extremite_aval": "0 0",
        },
        {
            "id": "Zone_F",
            "id_eltamont": "C_to_D",
            "id_eltaval": "Ext_Out",
            "extremite_amont": f"0 {2 * radious}",
            "extremite_aval": f"{-np. pi * radious} {2 * radious}",
        },
    )

    links = ET.SubElement(network, "TRONCONS")

    points_circle = [(0, 0), (radious, radious), (0, 2 * radious), (-radious, radious)]
    circle_center = np.add.reduce([np.array(v) for v in points_circle]) / 4

    internal_points = create_internal_points(points_circle, circle_center)

    TP_LINKS_INTPTS = ("Zone_B", "Zone_C", "Zone_D", "Zone_E")
    DCT_LINKS_INTPTS = dict(zip(TP_LINKS_INTPTS, internal_points))

    for tr in TP_LINKS:
        link = ET.SubElement(links, "TRONCON")
        link = fix_values(link, tr)
        if tr.get("id") in DCT_LINKS_INTPTS.keys():
            intpoints = ET.SubElement(link, "POINTS_INTERNES")
            for x, y in DCT_LINKS_INTPTS.get(tr.get("id")):
                ET.SubElement(intpoints, "POINT_INTERNE", coordonnees=f"{x} {y}")

    # Adding Traffic assingment
    assignments = ET.SubElement(network, "CONNEXIONS")
    distributors = ET.SubElement(assignments, "REPARTITEURS")

    for dist in TP_CONNECTIONS:
        distributor = ET.SubElement(distributors, "REPARTITEUR")
        distributor.set("id", dist)
        #TODO: Add traffic assignment for each distributor 
        #TODO: Create network graph directed graph first and then infer from there
    

    endings = ET.SubElement(network, "EXTREMITES")
    for ending in TP_EXTREMES:
        ET.SubElement(endings, "EXTREMITE", id=ending)

    # Adding Single Route
    routes = ET.SubElement(network, "ROUTES")
    route = ET.SubElement(routes, "ROUTE", description="In to Out", id="R01")
    route_path = ET.SubElement(route, "TRONCONS")

    laps = 5
    TP_ROUTE = ("Zone_A",) + ("Zone_B", "Zone_C", "Zone_D", "Zone_E") * laps + ("Zone_F",)

    for lnk in TP_ROUTE:
        ET.SubElement(route_path, "TRONCON", id=lnk)

    tree = ET.ElementTree(root)
    return tree


if __name__ == "__main__":
    print("Creating XML internal points")

    create_xml_simulation()
    # points = extract_points("ring_200.xml")
    # #     center = points[0][0], points[1][1] # Semi circle
    # center = np.add.reduce(points) / 4
    # print(f"center at: {center}")
    # int_points = create_internal_points(points, center, n_points=18, k=4)
    # xml_files = []
    # for i, int_point_seg in enumerate(int_points):
    #     xml_files = create_xml_data(int_point_seg)
    #     xml_files.write(f"filename_{i}.xml", pretty_print=True)
