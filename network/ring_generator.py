import os
import numpy as np
from lxml import etree as ET
import networkx as nx
import sys


def RotMat(theta: float) -> np.array:
    """Produces a rotation matrix 
    
    :param theta: angle to rotate
    :type theta: float
    :return: rotation matrix
    :rtype: np.array
    """
    return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])


def RotateVector(p: np.array, theta: float) -> np.array:
    """Rotates a vector p an angle theta counterclockwise
    
    :param p: vector to rotate
    :type p: np.array
    :param theta: angle in radians
    :type theta: float
    :return: vector once rotation is applied
    :rtype: np.array
    """
    return RotMat(theta) @ p


def RotationsVector(p: np.array, n_max: int = 15, k: int = 4) -> list:
    """Computes all vectors corresponding in steps of 2\pi / (k n_max)
    
    :param p: vector
    :type p: np.array
    :param n_max: max number of partitions, defaults to 10
    :type n_max: int, optional
    :param k: upper bound to rotate, defaults to 4 or pi/2
    :type k: int, optional
    :return: list of all vectors with corresponding rotations 
    :rtype: list
    """
    l_rot = []
    for n in range(0, n_max):
        theta = 2 * np.pi * n / (n_max * k)
        l_rot.append(RotateVector(p, theta))
    return l_rot


def create_internal_points(points: list, center: tuple, n_points: int = 15, k=4) -> tuple:
    """Compute a set of internal points for a ring-road 
    
    
    :return: list of points within the circle
    :rtype: tuple 
    """

    vct_pts = [np.array(pt) for pt in points]
    center = np.array(center)
    vct_ctr = [pt - center for pt in vct_pts]
    int_points = [RotationsVector(v, n_points, k)[1:] + center for v in vct_ctr]

    return int_points


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


def fix_values(element: ET.ElementTree, dct_values: dict) -> ET.ElementTree:
    """ Fix a set of values within a dictionary as attributes of an XML subelement
    
    :param element: XML tree element
    :type element: ET.ElementTree
    :param dct_values: attributes to fix
    :type dct_values: dict
    :return: Tranformed XML tree element
    :rtype: ET.ElementTree
    """
    for k, v in dct_values.items():
        element.set(k, v)
    return element


def create_line_graph(graph: nx.DiGraph) -> nx.DiGraph:
    """ Exchanges nodes for edges, and edges for links. 
        Returns the line graph of the graph or digraph. 
        Customized version of networkx.line_graph
    
    :param graph: Directed Graph 
    :type graph: nx.DiGraph
    :return: Directed graph where
    :rtype: nx.DiGraph
    """
    L = nx.DiGraph()
    for from_node in graph.edges():
        fromnode_name = graph.get_edge_data(*from_node).get("id")
        for to_node in graph.out_edges(from_node[1]):
            tonode_name = graph.get_edge_data(*to_node).get("id")
            L.add_edge(fromnode_name, tonode_name, id=from_node[1])
    return L


def transform_dict_data(orig_dct: dict, key: str):
    """ Type dispatcher to modify get method in terms of return type 

    :param orig_dct: original dictionary 
    :type orig_dct: dict
    :param key: string with key to transform
    :type key: str
    :return: value
    :rtype: matches KEYTYPEMAP
    """
    nested = lambda x: abs(float(x))
    KEYTYPEMAP = {"id": str, "w": nested, "kx": float, "vx": float, "N": int, "C": float}
    return KEYTYPEMAP.get(key, orig_dct.get(key))(orig_dct.get(key))


def create_xml_simulation(N_HDV: int, N_CAV: int = 2) -> ET.ElementTree:

    # --------------------------------------------------------------
    # Traffic simulation parameters

    # Scenario design
    # Sequence of vehicles
    SEQ_1 = (("HDV", N_HDV // 2), ("CAV", 1), ("HDV", N_HDV // 2), ("CAV", N_CAV - 1), ("CAV", 0))

    # Simulation info
    DCT_SIMULATION = {
        "id": "simID",
        "pasdetemps": "0.3",
        "debut": "07:00:00",
        "fin": "07:15:00",
        "loipoursuite": "exacte",
        "comportementflux": "iti",
        "date": "2019-06-19",
        "titre": "ring_simulation",
        "proc_deceleration": "false",
        "seed": "1",
    }
    DCT_EXPORT = {
        "trace_route": "false",
        "trajectoires": "true",
        "debug": "false",
        "debug_matrice_OD": "false",
        "debug_SAS": "false",
        "csv": "true",
    }
    DCT_TRAFIC = {"id": "trafID", "accbornee": "true", "coeffrelax": "4", "chgtvoie_ghost": "false"}
    DCT_NETWORK_INFO = {"id": "resID"}
    DCT_SCENARIO = {
        "id": "defaultScenario",
        "simulation_id": DCT_SIMULATION.get("id"),
        "trafic_id": DCT_TRAFIC.get("id"),
        "reseau_id": DCT_NETWORK_INFO.get("id"),
        "dirout": "data",
        "prefout": "simout_ring_data",
    }

    # Vehicle data
    TP_VEHTYPES = (
        {"id": "HDV", "w": "-5", "kx": "0.12", "vx": "25"},
        {"id": "CAV", "w": "-5", "kx": "0.12", "vx": "25"},
    )
    TP_ACCEL = (
        {"ax": "1.5", "vit_sup": "5.8"},
        {"ax": "1", "vit_sup": "8"},
        {"ax": "0.5", "vit_sup": "infini"},
    )

    # Road network data
    TP_TRONCONS = ("Zone_A", "Zone_B", "Zone_C", "Zone_D", "Zone_E", "Zone_F")
    TP_EXTREMES = ("Ext_In", "Ext_Out")
    TP_CONNECTIONS = ("A_to_B", "B_to_C", "C_to_D", "D_to_E")
    TP_LINKS_INTPTS = ("Zone_B", "Zone_C", "Zone_D", "Zone_E")

    radious = 200

    TP_LINKS = (
        {
            "id": "Zone_A",
            "id_eltamont": "Ext_In",
            "id_eltaval": "A_to_B",
            "extremite_amont": f"-{np.pi * radious} 0",
            "extremite_aval": "0 0",
            "largeur_voie": "3",
        },
        {
            "id": "Zone_B",
            "id_eltamont": "A_to_B",
            "id_eltaval": "B_to_C",
            "extremite_amont": "0 0",
            "extremite_aval": f"{radious} {radious}",
            "largeur_voie": "3",
        },
        {
            "id": "Zone_C",
            "id_eltamont": "B_to_C",
            "id_eltaval": "C_to_D",
            "extremite_amont": f"{radious} {radious}",
            "extremite_aval": f"0 {2 * radious}",
            "largeur_voie": "3",
        },
        {
            "id": "Zone_D",
            "id_eltamont": "C_to_D",
            "id_eltaval": "D_to_E",
            "extremite_amont": f"0 {2 * radious}",
            "extremite_aval": f"{-radious} {radious}",
            "largeur_voie": "3",
        },
        {
            "id": "Zone_E",
            "id_eltamont": "D_to_E",
            "id_eltaval": "A_to_B",
            "extremite_amont": f"{-radious} {radious}",
            "extremite_aval": "0 0",
            "largeur_voie": "3",
        },
        {
            "id": "Zone_F",
            "id_eltamont": "C_to_D",
            "id_eltaval": "Ext_Out",
            "extremite_amont": f"0 {2 * radious}",
            "extremite_aval": f"{-np. pi * radious} {2 * radious}",
            "largeur_voie": "3",
        },
    )

    # Creation of Routing:
    TP_CROSSINGS = TP_EXTREMES + TP_CONNECTIONS
    TP_LINKAGE = tuple(
        (
            road.get("id_eltamont"),
            road.get("id_eltaval"),
            {
                "id": road.get("id"),
                "extremite_amont": road.get("extremite_amont"),
                "extremite_aval": road.get("extremite_aval"),
                "largeur_voie": "3",
            },
        )
        for road in TP_LINKS
    )

    # Number of laps
    laps = 5
    TP_ROUTE = (
        ("Zone_A",)
        + ("Zone_B", "Zone_C", "Zone_D", "Zone_E") * laps
        + ("Zone_B", "Zone_C", "Zone_F")
    )

    # Demand Data
    TP_N = ({"N": N_CAV}, {"N": N_HDV})
    td = transform_dict_data
    DCT_VEHTYPES_PAR = {
        td(veh, "id"): {
            "kx": td(veh, "kx"),
            "vx": td(veh, "vx"),
            "w": td(veh, "w"),
            "C": (td(veh, "kx") * td(veh, "vx") * td(veh, "w")) / (td(veh, "vx") + td(veh, "w")),
            "N": td(nveh, "N"),
        }
        for veh, nveh in zip(TP_VEHTYPES, TP_N)
    }
    DCT_NW_DATA = {k: {"kc": d["C"] / d["vx"]} for k, d in DCT_VEHTYPES_PAR.items()}
    # Update wo assignment
    [y.update(DCT_NW_DATA.get(x)) for x, y in DCT_VEHTYPES_PAR.items()]

    # Scenario design
    rep_fnc = lambda d: "0 1" if d == "CAV" else "1 0"
    TP_DEMAND_TIME = list(
        [
            DCT_VEHTYPES_PAR.get(key).get("C"),
            value / DCT_VEHTYPES_PAR.get(key).get("C"),
            rep_fnc(key),
        ]
        for key, value in SEQ_1
    )

    # Infrastructure design
    radious = np.add.reduce(
        [veh.get("N") / (veh.get("kc") * np.pi * 2) for veh in DCT_VEHTYPES_PAR.values()]
    )
    print(f"R: {radious}")
    print(f"Circle length: {2 * np.pi * radious}")

    # Tricking demand
    TP_DEMAND_TIME[-1][0] = 0
    TP_DEMAND_TIME[-1][1] = laps * 2 * np.pi * radious / 25  # estimated time 3 laps

    # --------------------------------------------------------------

    # Create a general scheme
    attr_qname = ET.QName("http://www.w3.org/2001/XMLSchema-instance", "noNamespaceSchemaLocation")
    nsmap = {"xsi": "http://www.w3.org/2001/XMLSchema-instance"}
    root = ET.Element("ROOT_SYMUBRUIT", {attr_qname: "reseau.xsd"}, version="2.05", nsmap=nsmap)

    # Fill simulation
    simulations = ET.SubElement(root, "SIMULATIONS")
    simulation = ET.SubElement(simulations, "SIMULATION")
    simulation = fix_values(simulation, DCT_SIMULATION)
    export = ET.SubElement(simulation, "RESTITUTION")
    export = fix_values(export, DCT_EXPORT)

    # Add traffic
    trafics = ET.SubElement(root, "TRAFICS")
    trafic = ET.SubElement(trafics, "TRAFIC")
    trafic = fix_values(trafic, DCT_TRAFIC)

    # Add troncon
    troncons = ET.SubElement(trafic, "TRONCONS")
    for tron in TP_TRONCONS:
        troncon = ET.SubElement(troncons, "TRONCON")
        troncon.set("id", tron)

    # Add vehicle_type
    vehtypes = ET.SubElement(trafic, "TYPES_DE_VEHICULE")
    for veh in TP_VEHTYPES:
        veh_type = ET.SubElement(vehtypes, "TYPE_DE_VEHICULE")
        veh_type = fix_values(veh_type, veh)
        acc_maps = ET.SubElement(veh_type, "ACCELERATION_PLAGES")
        # Add acceleration
        for acc in TP_ACCEL:
            acc_map = ET.SubElement(acc_maps, "ACCELERATION_PLAGE")
            acc_map = fix_values(acc_map, acc)

    # Adding borders + demand
    extremes = ET.SubElement(trafic, "EXTREMITES")

    for ext in TP_EXTREMES:
        extreme = ET.SubElement(extremes, "EXTREMITE")
        extreme.set("id", ext)
        if ext == "Ext_In":
            flux_global = ET.SubElement(extreme, "FLUX_GLOBAL")
            flux = ET.SubElement(flux_global, "FLUX")
            demands = ET.SubElement(flux, "DEMANDES")
            for level, time, _ in TP_DEMAND_TIME:
                demand = ET.SubElement(demands, "DEMANDE")
                demand.set("niveau", str(level))
                demand.set("duree", str(time))
            split_ratios = ET.SubElement(flux, "REP_DESTINATIONS")
            split_ratio = ET.SubElement(split_ratios, "REP_DESTINATION")
            destination = ET.SubElement(split_ratio, "DESTINATION")
            destination.set("coeffOD", "1")
            destination.set("sortie", "Ext_Out")
            route_dstn = ET.SubElement(destination, "ROUTE")
            route_dstn.set("coeffAffectation", "1")
            route_dstn.set("id", "R01")
            split_vehtypes = ET.SubElement(flux_global, "REP_TYPEVEHICULES")
            for _, time, trfassg in TP_DEMAND_TIME:
                split_vehtype = ET.SubElement(split_vehtypes, "REP_TYPEVEHICULE")
                split_vehtype.set("duree", str(time))
                split_vehtype.set("coeffs", trfassg)

    # Adding internal connections
    connections = ET.SubElement(trafic, "CONNEXIONS_INTERNES")
    for con in TP_CONNECTIONS:
        connection = ET.SubElement(connections, "CONNEXION_INTERNE")
        connection.set("id", con)

    # Adding network links + internal points
    networks = ET.SubElement(root, "RESEAUX")
    network = ET.SubElement(networks, "RESEAU")

    network = fix_values(network, DCT_NETWORK_INFO)

    links = ET.SubElement(network, "TRONCONS")

    # This is just an update based on the updated radious
    TP_LINKS = (
        {
            "id": "Zone_A",
            "id_eltamont": "Ext_In",
            "id_eltaval": "A_to_B",
            "extremite_amont": f"-{np.pi * radious} 0",
            "extremite_aval": "0 0",
            "largeur_voie": "3",
        },
        {
            "id": "Zone_B",
            "id_eltamont": "A_to_B",
            "id_eltaval": "B_to_C",
            "extremite_amont": "0 0",
            "extremite_aval": f"{radious} {radious}",
            "largeur_voie": "3",
        },
        {
            "id": "Zone_C",
            "id_eltamont": "B_to_C",
            "id_eltaval": "C_to_D",
            "extremite_amont": f"{radious} {radious}",
            "extremite_aval": f"0 {2 * radious}",
            "largeur_voie": "3",
        },
        {
            "id": "Zone_D",
            "id_eltamont": "C_to_D",
            "id_eltaval": "D_to_E",
            "extremite_amont": f"0 {2 * radious}",
            "extremite_aval": f"{-radious} {radious}",
            "largeur_voie": "3",
        },
        {
            "id": "Zone_E",
            "id_eltamont": "D_to_E",
            "id_eltaval": "A_to_B",
            "extremite_amont": f"{-radious} {radious}",
            "extremite_aval": "0 0",
            "largeur_voie": "3",
        },
        {
            "id": "Zone_F",
            "id_eltamont": "C_to_D",
            "id_eltaval": "Ext_Out",
            "extremite_amont": f"0 {2 * radious}",
            "extremite_aval": f"{-np. pi * radious} {2 * radious}",
            "largeur_voie": "3",
        },
    )

    points_circle = [(0, 0), (radious, radious), (0, 2 * radious), (-radious, radious)]
    circle_center = np.add.reduce([np.array(v) for v in points_circle]) / 4

    internal_points = create_internal_points(points_circle, circle_center)

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

    # Network Creation
    RoadGraph = nx.DiGraph()
    RoadGraph.add_nodes_from(TP_CROSSINGS)
    RoadGraph.add_edges_from(TP_LINKAGE)

    # Exchange in between RoadNetwork and Linkage network
    LinkageGraph = create_line_graph(RoadGraph)

    for edge in LinkageGraph.edges():
        distributor = ET.SubElement(distributors, "REPARTITEUR")
        node = LinkageGraph.get_edge_data(*edge).get("id")
        distributor.set("id", node)
        ath_mvts = ET.SubElement(distributor, "MOUVEMENTS_AUTORISES")
        ath_mvt = ET.SubElement(ath_mvts, "MOUVEMENT_AUTORISE")
        ath_mvt.set("id_troncon_amont", edge[0])
        out_mvts = ET.SubElement(ath_mvt, "MOUVEMENT_SORTIES")
        out_mvt = ET.SubElement(out_mvts, "MOUVEMENT_SORTIE")
        out_mvt.set("id_troncon_aval", edge[1])

    endings = ET.SubElement(assignments, "EXTREMITES")
    for ending in TP_EXTREMES:
        ET.SubElement(endings, "EXTREMITE", id=ending)

    # Adding Single Route
    routes = ET.SubElement(network, "ROUTES")
    route = ET.SubElement(routes, "ROUTE", description="In to Out", id="R01")
    route_path = ET.SubElement(route, "TRONCONS_")

    for lnk in TP_ROUTE:
        ET.SubElement(route_path, "TRONCON_", id=lnk)

    # Adding Scenario
    scenarios = ET.SubElement(root, "SCENARIOS")
    scenario = ET.SubElement(scenarios, "SCENARIO")
    scenario = fix_values(scenario, DCT_SCENARIO)

    tree = ET.ElementTree(root)
    return tree


if __name__ == "__main__":
    print("Creating XML internal points")

    N = int(sys.argv[1])
    x1 = create_xml_simulation(N)
    x1.write(f"ring_{N}_cav.xml", pretty_print=True, encoding="UTF8")

