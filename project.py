import graph_tool.all as gt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from geopy.distance import great_circle
import numpy as np
import argparse
import os
import heapq
import csv
from pathlib import Path

DELAY_PER_KM = 4.8

# алгоритм Дейкстры для поиска максимального пути при данном расположении контроллера
# улучшения не использовал, так как они, на мой взгляд, не могут оптимизировать выполнение данной задачи
def Dijkstra(v, g):
    times = [float("inf")] * g.num_vertices()
    times[int(v)] = 0
    visited = set()
    queue = [(0, v)]
    while queue:
        u_time, u = heapq.heappop(queue)
        if int(u) in visited:
            continue
        visited.add(int(u))
        for e in u.all_edges():
            w = e.target() if u == e.source() else e.source()
            new_w_time = u_time + g.ep.distance[e]
            if new_w_time < times[int(w)]:
                times[int(w)] = new_w_time
                heapq.heappush(queue, (new_w_time, w))
    return np.max(times)

# алгоритм Прима для поиска минимального остовного дерева
def Prim(v, g):
    visited = {int(v)}
    for _ in range(g.num_vertices() - 1):
        min_weight = float("inf")
        for u in visited:
            for e in g.vertex(u).all_edges():
                w = e.target() if u == int(e.source()) else e.source()
                if g.ep.distance[e] < min_weight and int(w) not in visited:
                    chosen_v = w
                    chosen_e = e
                    min_weight = g.ep.distance[e]
        visited.add(int(chosen_v))
        g.ep.color[chosen_e] = "blue"
    return g

# функция поиска пути до вершин в дереве
def FindPath(v, visited, paths, delays, g):
    for u in v.all_neighbors():
        if int(u) in visited:
            continue
        paths[int(u)] = paths[int(v)] + [int(u) + 1]
        e = g.edge(u,v)
        delays[int(u)] = delays[int(v)] + g.ep.distance[e] * DELAY_PER_KM
        visited.add(int(u))
        FindPath(u, visited, paths, delays, g)
    return

# получение аргументов
parser = argparse.ArgumentParser()

parser.add_argument(
    "-t",
    help = "Файл с топологией сети",
    metavar = '',
    required = True
)
parser.add_argument(
    "-k",
    type = int,
    help = "Критерий выбора (1 или 2)",
    metavar = '',
    choices = [1, 2],
    default = 1
)
parser.add_argument(
    "-p",
    help = "Директория, где будут размещены результаты работы",
    metavar = '',
    default = 'tests'
)
args = parser.parse_args()

fl = args.t
if not os.path.exists(fl):
    raise TypeError(f"Файл '{fl}' не существует")
if not fl.lower().endswith((".gml", ".graphml", ".xml")):
    raise TypeError("Формат файла должен быть .gml или .graphml")

path = Path(args.p)
if not path.exists():
    raise TypeError(f"Директории '{path}' не существует")
if not path.is_dir():
    raise TypeError(f"'{path}' не является директорией")
if not os.access(str(path), os.W_OK):
    raise TypeError(f"Нет доступа к директории '{path}'")
path = args.p

g = gt.load_graph(fl)
g.set_directed(False)

# проверим наличие необходимых свойств у вершин
atributes = {"label", "latitude", "longitude"}
for atribute in atributes:
    if atribute not in g.vertex_properties:
        raise ValueError(f"Необходим атрибут '{atribute}' для вершин графа")

# проверим связность графа
labels, _ = gt.label_components(g)
if len(set(labels)) != 1:
    raise ValueError("Граф не связный")

lats = g.vertex_properties["latitude"]  # широта
lons = g.vertex_properties["longitude"] # долгота

# создаём свойства цвета для рёбер и вершин
g.ep["color"] = g.new_edge_property("string", val = "white")
g.vp["color"] = g.new_vertex_property("string", val = "yellow")

# создаём свойство веса (задержки) для рёбер
distance = g.new_edge_property("double")

for e in g.edges():
    src, tgt = e.source(), e.target()
    distance[e] = great_circle((lats[src], lons[src]), (lats[tgt], lons[tgt])).km
    # используется упрощённая функция подсчёта расстояния из библиотеки geopy
    # источник: https://geopy.readthedocs.io/en/stable/#module-geopy.distance

g.ep["distance"] = distance

if args.k == 1:
    min_time = float("inf")
    for v in g.vertices():
        v_time = Dijkstra(v, g)
        if v_time < min_time:
            min_time = v_time
            controller = v

    g.vp.color[controller] = "red"
else:
    # функция изменит цвет выбранных для дерева рёбер, после используем фильтр
    filt_g = gt.GraphView(Prim(g.vertex(0), g), efilt = lambda e: g.ep.color[e] == "blue")
    min_time = float("inf")
    for v in filt_g.vertices():
        v_time = Dijkstra(v, filt_g)
        if v_time < min_time:
            min_time = v_time
            controller = v

    g.vp.color[controller] = "red"

# создаём границы для отображения
padding = 2.0 
lat_min, lat_max = np.min(lats) - padding, np.max(lats) + padding
lon_min, lon_max = np.min(lons) - padding, np.max(lons) + padding



# создаём карту
fig = plt.figure(figsize=(15, 12))
ax = fig.add_subplot(1, 1, 1, projection = ccrs.PlateCarree())
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs = ccrs.PlateCarree())
ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor = 'lightgreen')
ax.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor = 'lightblue')
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth = 0.8)
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle = ':', linewidth = 0.5)
ax.gridlines(draw_labels = True, linestyle = '--', alpha = 0.7)


scale = max(lon_max - lon_min, lat_max - lat_min)

# ребра
for e in g.edges():
    src, tgt = e.source(), e.target()
    ax.plot(
        [lons[src], lons[tgt]],
        [lats[src], lats[tgt]],
        color = g.ep.color[e],
        linewidth = scale / 25,
        transform = ccrs.Geodetic()
    )

# вершины
for v in g.vertices():
    ax.plot(
        lons[v], lats[v], 
        color = g.vp.color[v],
        marker = 'o', 
        markersize = scale / 5,
        transform = ccrs.Geodetic()
    )

# вывод в csv файлы
fl_topo = path + '/' + Path(args.t).stem + "_topo.csv"
with open(fl_topo, 'w', newline='') as csv_topo:
    writer = csv.writer(csv_topo)
    writer.writerow([
        "Node 1\n(id)",
        "Node 1\n(label)",
        "Node 1\n(longitude)",
        "Node 1\n(latitude)",
        "Node 2\n(id)",
        "Node 2\n(label)",
        "Node 2\n(longitude)",
        "Node 2\n(latitude)",
        "Distance\n(km)",
        "Delay\n(mks)"
    ])
    
    for i in range(g.num_vertices()):
        v = g.vertex(i)
        for e in v.all_edges():
            u = e.target() if v == e.source() else e.source()
            writer.writerow([
                int(v) + 1,
                g.vp.label[v],
                g.vp.longitude[v],
                g.vp.latitude[v],
                int(u) + 1,
                g.vp.label[u],
                g.vp.longitude[u],
                g.vp.latitude[u],
                g.ep.distance[e],
                g.ep.distance[e] * DELAY_PER_KM
            ])

if args.k == 2:
    fl_routes = path + '/' + Path(args.t).stem + "_routes.csv"
    with open(fl_routes, 'w', newline='') as csv_routes:
        writer = csv.writer(csv_routes)
        writer.writerow([
            "Node 1 (id)",
            "Node 2 (id)",
            "Path type",
            "Path",
            "Delay (mks)"
        ])

        visited = {int(controller)}
        paths = [[] for _ in range(filt_g.num_vertices())]
        paths[int(controller)] = [int(controller) + 1]
        delays = [0] * filt_g.num_vertices()
        FindPath(controller, visited, paths, delays, filt_g)

        for v in filt_g.vertices():
            if (v == controller):
                continue
            writer.writerow([
                int(controller) + 1,
                int(v) + 1,
                "main",
                str(paths[int(v)]),
                delays[int(v)]
            ])

plt.show()
fig.savefig(path + '/' + "map.png", dpi = 300, bbox_inches = "tight", facecolor = "white")