"""
Valencia Flood Evacuation Routing — v2 (Road Network + Two-Phase)
==================================================================
QGIS Python Console veya Script Editörü'nde çalıştırın.

SENARYO:
  Faz 1 — Halk (yaya/araç) → Yol ağı üzerinden → Toplanma noktası
  Faz 2 — Devlet araçları  → Yol ağı üzerinden → Toplanma noktası → Güvenli tesis

Gereksinimler:
  - QGIS 3.x
  - 'valenciahrsl' ve 'valenciariskflood' katmanları yüklü olmalı (EPSG:25830)
  - İnternet bağlantısı (Overpass API)
  - 'requests' kütüphanesi (genellikle QGIS'te var; yoksa: pip install requests)
"""

import heapq
import math
import random
from collections import defaultdict, deque

import numpy as np
import requests
from osgeo import gdal

from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsFeature,
    QgsField,
    QgsFields,
    QgsGeometry,
    QgsPointXY,
    QgsProject,
    QgsCategorizedSymbolRenderer,
    QgsRendererCategory,
    QgsLineSymbol,
    QgsMarkerSymbol,
    QgsSingleSymbolRenderer,
    QgsVectorLayer,
    QgsWkbTypes,
)
from qgis.PyQt.QtCore import QVariant

# ══════════════════════════════════════════════════════════════
# 0. AYARLAR
# ══════════════════════════════════════════════════════════════
POPULATION_LAYER_NAME  = "valenciahrsl"
FLOODRISK_LAYER_NAME   = "valenciariskflood"

# Overlap eşiği: üst kaç % piksel kaynak kabul edilsin
SOURCE_PERCENTILE      = 95

# Flood risk eşiği — hedef/toplanma noktasının flood riski bu değerin altında olmalı
SAFE_RISK_THRESHOLD    = 0.10

# Dijkstra kenar maliyet çarpanı: yüksek → riskli yollardan daha fazla kaçınır
FLOOD_WEIGHT           = 20.0

# Minimum cluster boyutu (piksel)
MIN_CLUSTER_SIZE       = 5

# OSM sorgu/bekleme süresi (saniye)
OSM_TIMEOUT            = 90

OVERPASS_URL           = "https://overpass-api.de/api/interpreter"

# Her küme için kaç alternatif Assembly Point denensin?
ASSEMBLY_CANDIDATES    = 3

# Yol tipine göre hız düzeltme çarpanı (düşük → hızlı, yüksek → yavaş)
ROAD_SPEED = {
    "motorway": 0.5, "motorway_link": 0.6,
    "trunk": 0.6,    "trunk_link": 0.65,
    "primary": 0.7,  "primary_link": 0.75,
    "secondary": 0.8,"secondary_link": 0.85,
    "tertiary": 0.9, "tertiary_link": 0.95,
    "residential": 1.1, "living_street": 1.3,
    "service": 1.4,   "track": 2.0,
    "path": 2.5,      "footway": 2.5,
    "steps": 4.0,     "unclassified": 1.2,
}
DEFAULT_ROAD_SPEED = 1.5   # bilinmeyen tipler için


# ══════════════════════════════════════════════════════════════
# 1. RASTER YARDIMCILARI
# ══════════════════════════════════════════════════════════════

def get_layer(name):
    layers = QgsProject.instance().mapLayersByName(name)
    if not layers:
        raise ValueError(f"'{name}' katmanı bulunamadı.")
    return layers[0]


def raster_to_array(qlayer):
    ds = gdal.Open(qlayer.source())
    if ds is None:
        raise IOError(f"GDAL açamadı: {qlayer.source()}")
    band   = ds.GetRasterBand(1)
    arr    = band.ReadAsArray().astype(np.float32)
    nodata = band.GetNoDataValue()
    if nodata is not None:
        arr[arr == nodata] = 0.0
    gt      = ds.GetGeoTransform()
    crs_wkt = ds.GetProjection()
    ds = None
    return arr, gt, crs_wkt


def normalize(arr):
    mn, mx = arr.min(), arr.max()
    if mx == mn:
        return np.zeros_like(arr)
    return (arr - mn) / (mx - mn)


def pixel_to_geo(gt, row, col):
    x = gt[0] + col * gt[1] + row * gt[2]
    y = gt[3] + col * gt[4] + row * gt[5]
    return x, y


def geo_to_pixel(gt, x, y):
    """Geo koordinatı piksel (row, col)'a çevirir, sınır dışıysa None döner."""
    det = gt[1] * gt[5] - gt[2] * gt[4]
    if abs(det) < 1e-12:
        return None
    col = (gt[5] * (x - gt[0]) - gt[2] * (y - gt[3])) / det
    row = (gt[1] * (y - gt[3]) - gt[4] * (x - gt[0])) / det
    return int(row), int(col)


def sample_risk(risk_norm, gt, x, y):
    """Koordinatta flood risk değerini örnekler; sınır dışıysa 1.0 (çok riskli) döner."""
    shape = risk_norm.shape
    rc    = geo_to_pixel(gt, x, y)
    if rc is None:
        return 1.0
    r, c = rc
    if 0 <= r < shape[0] and 0 <= c < shape[1]:
        return float(risk_norm[r, c])
    return 1.0


# ══════════════════════════════════════════════════════════════
# 2. KAYNAK KÜMELERİ (yüksek nüfus × yüksek risk)
# ══════════════════════════════════════════════════════════════

def find_source_clusters(pop_arr, risk_arr, percentile=SOURCE_PERCENTILE):
    overlap  = normalize(pop_arr) * normalize(risk_arr)
    pos      = overlap[overlap > 0]
    if len(pos) == 0:
        return [], overlap
    threshold = np.percentile(pos, percentile)
    source_mask = overlap >= threshold

    rows, cols = source_mask.shape
    visited    = np.zeros_like(source_mask, dtype=bool)
    clusters   = []

    def bfs(sr, sc):
        cl  = []
        q   = deque([(sr, sc)])
        visited[sr, sc] = True
        while q:
            r, c = q.popleft()
            cl.append((r, c))
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if dr == 0 and dc == 0:
                        continue
                    nr, nc = r + dr, c + dc
                    if 0 <= nr < rows and 0 <= nc < cols:
                        if source_mask[nr, nc] and not visited[nr, nc]:
                            visited[nr, nc] = True
                            q.append((nr, nc))
        return cl

    for r in range(rows):
        for c in range(cols):
            if source_mask[r, c] and not visited[r, c]:
                cl = bfs(r, c)
                if len(cl) >= MIN_CLUSTER_SIZE:
                    clusters.append(cl)

    print(f"[Cluster] {len(clusters)} kaynak bölge bulundu "
          f"(overlap ≥ {threshold:.4f}, eşik %{percentile})")
    return clusters, overlap


def cluster_centroid_geo(cluster, gt):
    """Kümenin medyan pikselinin coğrafi koordinatını döndürür."""
    rows = [p[0] for p in cluster]
    cols = [p[1] for p in cluster]
    mr   = int(np.median(rows))
    mc   = int(np.median(cols))
    return pixel_to_geo(gt, mr, mc)


# ══════════════════════════════════════════════════════════════
# 3. CRS DÖNÜŞÜMÜ (EPSG:25830 ↔ WGS84)
# ══════════════════════════════════════════════════════════════

def make_transforms():
    src = QgsCoordinateReferenceSystem("EPSG:25830")
    wgs = QgsCoordinateReferenceSystem("EPSG:4326")
    proj = QgsProject.instance()
    to_wgs   = QgsCoordinateTransform(src, wgs, proj)
    from_wgs = QgsCoordinateTransform(wgs, src, proj)
    return to_wgs, from_wgs


def bbox_to_wgs84(bbox_25830, to_wgs):
    xmin, ymin, xmax, ymax = bbox_25830
    sw = to_wgs.transform(QgsPointXY(xmin, ymin))
    ne = to_wgs.transform(QgsPointXY(xmax, ymax))
    return sw.y(), sw.x(), ne.y(), ne.x()   # S, W, N, E


# ══════════════════════════════════════════════════════════════
# 4. OSM YOL AĞINI ÇEKME VE GRAF OLUŞTURMA
# ══════════════════════════════════════════════════════════════

def fetch_osm_roads(bbox_wgs84):
    """
    Overpass'tan yol ağını çeker.
    Döndürür: nodes dict {id: (lon, lat)}, ways list of {nodes:[id,...], tags:{}}
    """
    S, W, N, E = bbox_wgs84
    bbox_str   = f"{S},{W},{N},{E}"

    # Yalnızca araç/yaya geçişine uygun yollar
    query = f"""
[out:json][timeout:{OSM_TIMEOUT}];
(
  way["highway"~"^(motorway|trunk|primary|secondary|tertiary|
    unclassified|residential|living_street|service|
    motorway_link|trunk_link|primary_link|secondary_link|tertiary_link|
    pedestrian|footway|path|track)$"]
    ({bbox_str});
);
out body;
>;
out skel qt;
"""
    print(f"[OSM/Roads] Yol ağı sorgulanıyor ({bbox_str})...")
    try:
        r = requests.post(OVERPASS_URL, data={"data": query}, timeout=OSM_TIMEOUT + 15)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"[OSM/Roads] HATA: {e}")
        return {}, []

    nodes = {}
    ways  = []
    for el in data.get("elements", []):
        if el["type"] == "node":
            nodes[el["id"]] = (el["lon"], el["lat"])
        elif el["type"] == "way":
            ways.append({"nodes": el.get("nodes", []),
                         "tags":  el.get("tags", {})})

    print(f"[OSM/Roads] {len(nodes)} düğüm, {len(ways)} yol segmenti alındı.")
    return nodes, ways


def build_road_graph(osm_nodes, osm_ways, from_wgs, risk_norm, gt):
    """
    OSM düğüm/yollarından ağırlıklı çift yönlü graf oluşturur.
    Düğüm ID'si: OSM node id (int)
    Düğüm konumu: EPSG:25830 (x, y)
    Kenar ağırlığı: mesafe × yol_hızı × (1 + FLOOD_WEIGHT × ortalama_risk)

    Döndürür:
      graph  : {node_id: [(neighbour_id, cost), ...]}
      pos    : {node_id: (x25830, y25830)}
    """
    # WGS84 → EPSG:25830
    pos = {}
    for nid, (lon, lat) in osm_nodes.items():
        pt = from_wgs.transform(QgsPointXY(lon, lat))
        pos[nid] = (pt.x(), pt.y())

    graph = defaultdict(list)

    for way in osm_ways:
        nids     = way["nodes"]
        tags     = way["tags"]
        hw_type  = tags.get("highway", "unclassified")
        one_way  = tags.get("oneway", "no") in ("yes", "1", "true")
        spd_mult = ROAD_SPEED.get(hw_type, DEFAULT_ROAD_SPEED)

        for i in range(len(nids) - 1):
            a, b = nids[i], nids[i + 1]
            if a not in pos or b not in pos:
                continue
            ax, ay = pos[a]
            bx, by = pos[b]
            dist   = math.hypot(bx - ax, by - ay)   # metre (EPSG:25830)
            if dist < 0.1:
                continue

            # Ortalama flood riski bu segmentte
            risk_a = sample_risk(risk_norm, gt, ax, ay)
            risk_b = sample_risk(risk_norm, gt, bx, by)
            avg_r  = (risk_a + risk_b) / 2.0

            cost   = dist * spd_mult * (1.0 + FLOOD_WEIGHT * avg_r)

            graph[a].append((b, cost))
            if not one_way:
                graph[b].append((a, cost))

    print(f"[Graph] Yol grafı oluşturuldu: {len(graph)} düğüm, "
          f"{sum(len(v) for v in graph.values())} kenar")
    return dict(graph), pos


def nearest_node(px, py, pos, max_dist=2000):
    """
    (px, py) koordinatına en yakın graf düğümünü döndürür.
    max_dist metre ötesindeyse None döner.
    """
    best_id   = None
    best_dist = float("inf")
    for nid, (x, y) in pos.items():
        d = math.hypot(x - px, y - py)
        if d < best_dist:
            best_dist = d
            best_id   = nid
    if best_dist > max_dist:
        return None, best_dist
    return best_id, best_dist


# ══════════════════════════════════════════════════════════════
# 5. OSM TESİS SORGULAMA
# ══════════════════════════════════════════════════════════════

def fetch_osm_facilities(bbox_wgs84, from_wgs):
    """
    Toplanma noktaları ve güvenli tesisler:
      Faz 1 Hedef  → assembly_point, open_space, park, school, stadium
      Faz 2 Hedef  → hospital, fire_station, police, clinic, shelter
    Döndürür:
      assembly_pts : [(x25830, y25830, name, kind), ...]
      safe_pts     : [(x25830, y25830, name, kind), ...]
    """
    S, W, N, E = bbox_wgs84
    bbox_str   = f"{S},{W},{N},{E}"

    # Akademik kaynaklarda önerilen toplanma / tahliye tesisleri
    query = f"""
[out:json][timeout:{OSM_TIMEOUT}];
(
  node["emergency"="assembly_point"]({bbox_str});
  node["amenity"~"school|university|stadium|parking"]({bbox_str});
  node["leisure"~"park|sports_centre|stadium"]({bbox_str});
  node["amenity"~"hospital|fire_station|police|clinic|shelter"]({bbox_str});
  way["emergency"="assembly_point"]({bbox_str});
  way["amenity"~"school|university|stadium|parking"]({bbox_str});
  way["leisure"~"park|sports_centre|stadium"]({bbox_str});
  way["amenity"~"hospital|fire_station|police|clinic|shelter"]({bbox_str});
);
out center;
"""
    print(f"[OSM/Facilities] Tesisler sorgulanıyor...")
    try:
        r = requests.post(OVERPASS_URL, data={"data": query}, timeout=OSM_TIMEOUT + 15)
        r.raise_for_status()
        data = r.json()
    except Exception as e:
        print(f"[OSM/Facilities] HATA: {e}")
        return [], []

    ASSEMBLY_KINDS = {"assembly_point", "school", "university",
                      "stadium", "park", "sports_centre", "parking"}
    SAFE_KINDS     = {"hospital", "fire_station", "police",
                      "clinic", "shelter"}

    assembly_pts = []
    safe_pts     = []

    for el in data.get("elements", []):
        if el["type"] == "node":
            lon, lat = el["lon"], el["lat"]
        elif el["type"] == "way" and "center" in el:
            lon, lat = el["center"]["lon"], el["center"]["lat"]
        else:
            continue

        tags = el.get("tags", {})
        kind = (tags.get("emergency")
                or tags.get("amenity")
                or tags.get("leisure") or "unknown")
        name = tags.get("name", kind)

        pt   = from_wgs.transform(QgsPointXY(lon, lat))
        x25, y25 = pt.x(), pt.y()

        if kind in SAFE_KINDS:
            safe_pts.append((x25, y25, name, kind))
        elif kind in ASSEMBLY_KINDS:
            assembly_pts.append((x25, y25, name, kind))
        # Bazı tesisler her iki listeye de eklenebilir (örn. büyük otopark)

    print(f"[OSM/Facilities] {len(assembly_pts)} toplanma noktası, "
          f"{len(safe_pts)} güvenli tesis bulundu.")
    return assembly_pts, safe_pts


# ══════════════════════════════════════════════════════════════
# 6. YOLLARDA DİJKSTRA
# ══════════════════════════════════════════════════════════════

def dijkstra_graph(graph, sources_nodes, dest_nodes):
    """
    Yol grafı üzerinde Dijkstra.
    sources_nodes : [node_id, ...]
    dest_nodes    : set of node_id

    Döndürür: (path of node_ids, hit_dest_id) ya da (None, None)
    """
    dist = {}
    prev = {}
    heap = []

    for s in sources_nodes:
        if s in graph:
            dist[s] = 0.0
            heapq.heappush(heap, (0.0, s))

    while heap:
        d, u = heapq.heappop(heap)

        if u in dest_nodes:
            path = []
            cur  = u
            while cur in prev:
                path.append(cur)
                cur = prev[cur]
            path.append(cur)
            path.reverse()
            return path, u

        if d > dist.get(u, float("inf")):
            continue

        for (v, cost) in graph.get(u, []):
            nd = d + cost
            if nd < dist.get(v, float("inf")):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(heap, (nd, v))

    return None, None


# ══════════════════════════════════════════════════════════════
# 7. GÜVENLI NOKTA FİLTRELEME (flood risk kontrolü)
# ══════════════════════════════════════════════════════════════

def filter_safe_nodes(pts_list, pos, risk_norm, gt,
                      threshold=SAFE_RISK_THRESHOLD, max_snap_dist=2000):
    """
    OSM tesislerini yol düğümlerine yaklaştırır ve flood riskini filtreler.
    Döndürür: [(node_id, label, x, y), ...]
    """
    result = []
    for (x, y, name, kind) in pts_list:
        risk_val = sample_risk(risk_norm, gt, x, y)
        if risk_val >= threshold:
            continue   # Yüksek riskli tesis → geçersiz
        nid, dist = nearest_node(x, y, pos, max_snap_dist)
        if nid is None:
            continue
        # Tesise yakın düğümün de riski kontrol et
        nx, ny  = pos[nid]
        risk_nd = sample_risk(risk_norm, gt, nx, ny)
        if risk_nd >= threshold:
            continue
        result.append((nid, f"{kind}: {name}", x, y))

    return result


# ══════════════════════════════════════════════════════════════
# 8. QGIS KATMAN OLUŞTURMA
# ══════════════════════════════════════════════════════════════

def path_nodes_to_coords(path, pos):
    return [QgsPointXY(*pos[nid]) for nid in path if nid in pos]


def make_route_layer(routes, crs_str, layer_name, color_field="id"):
    """
    routes: list of dicts:
      {id, path_coords: [QgsPointXY], label, phase}
    """
    fields = QgsFields()
    fields.append(QgsField("id",    QVariant.Int))
    fields.append(QgsField("faz",   QVariant.Int))
    fields.append(QgsField("hedef", QVariant.String))

    vlayer = QgsVectorLayer(f"LineString?crs={crs_str}", layer_name, "memory")
    vlayer.startEditing()
    vlayer.dataProvider().addAttributes(fields)
    vlayer.updateFields()

    for route in routes:
        pts = route["path_coords"]
        if len(pts) < 2:
            continue
        feat = QgsFeature()
        feat.setGeometry(QgsGeometry.fromPolylineXY(pts))
        feat.setAttributes([route["id"], route["phase"], route["label"]])
        vlayer.addFeature(feat)

    vlayer.commitChanges()
    return vlayer


def make_point_layer(points, crs_str, layer_name):
    """
    points: [(x, y, label), ...]
    """
    fields = QgsFields()
    fields.append(QgsField("label", QVariant.String))

    vlayer = QgsVectorLayer(f"Point?crs={crs_str}", layer_name, "memory")
    vlayer.startEditing()
    vlayer.dataProvider().addAttributes(fields)
    vlayer.updateFields()

    for (x, y, label) in points:
        feat = QgsFeature()
        feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
        feat.setAttributes([label])
        vlayer.addFeature(feat)

    vlayer.commitChanges()
    return vlayer


def style_route_layer(vlayer, phase):
    """Faz 1 → turuncu, Faz 2 → kırmızı çizgiler"""
    color   = "255,140,0,255" if phase == 1 else "200,0,0,255"
    width   = "1.0"           if phase == 1 else "1.5"
    dash    = "5;3"           if phase == 1 else ""   # faz1 kesik
    sym     = QgsLineSymbol.createSimple({
        "color": color, "width": width,
        "capstyle": "round", "joinstyle": "round",
    })
    if dash:
        sym.symbolLayer(0).setCustomDashVector([5, 3])
        sym.symbolLayer(0).setUseCustomDashPattern(True)
    vlayer.setRenderer(QgsSingleSymbolRenderer(sym))
    vlayer.triggerRepaint()


def style_point_layer(vlayer, color_rgb, size="4"):
    sym = QgsMarkerSymbol.createSimple({
        "color": color_rgb, "size": size,
        "outline_color": "50,50,50,200",
        "outline_width": "0.3",
    })
    vlayer.setRenderer(QgsSingleSymbolRenderer(sym))
    vlayer.triggerRepaint()


# ══════════════════════════════════════════════════════════════
# 9. ANA ÇALIŞMA AKIŞI
# ══════════════════════════════════════════════════════════════

def run_evacuation_routing():
    print("=" * 65)
    print(" Valencia Flood Evacuation Routing  v2 — Road Network")
    print("=" * 65)

    # ── Raster okuma ──────────────────────────────────────────
    pop_layer  = get_layer(POPULATION_LAYER_NAME)
    risk_layer = get_layer(FLOODRISK_LAYER_NAME)

    pop_arr,  gt_pop,  _ = raster_to_array(pop_layer)
    risk_arr, gt_risk, _ = raster_to_array(risk_layer)

    # Boyut uyumu
    if pop_arr.shape != risk_arr.shape:
        from scipy.ndimage import zoom
        zf       = (pop_arr.shape[0] / risk_arr.shape[0],
                    pop_arr.shape[1] / risk_arr.shape[1])
        risk_arr = zoom(risk_arr, zf, order=1)
    gt = gt_pop

    risk_norm = normalize(risk_arr)
    print(f"[Raster] Pop boyutu: {pop_arr.shape}, Risk boyutu: {risk_arr.shape} ✓")

    # ── Kaynak kümeler ────────────────────────────────────────
    clusters, _ = find_source_clusters(pop_arr, risk_arr, SOURCE_PERCENTILE)
    if not clusters:
        print("[HATA] Kaynak bölge bulunamadı.")
        return

    # ── Bounding box ──────────────────────────────────────────
    ROWS, COLS = risk_arr.shape
    x0, y0 = pixel_to_geo(gt, 0, 0)
    x1, y1 = pixel_to_geo(gt, ROWS - 1, COLS - 1)
    bbox_25830 = (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))

    to_wgs, from_wgs = make_transforms()
    bbox_wgs84 = bbox_to_wgs84(bbox_25830, to_wgs)
    print(f"[Bbox] EPSG:25830 — {tuple(round(v,1) for v in bbox_25830)}")

    # ── OSM yol ağı ───────────────────────────────────────────
    osm_nodes, osm_ways = fetch_osm_roads(bbox_wgs84)
    if not osm_nodes:
        print("[HATA] OSM yol ağı alınamadı.")
        return

    graph, pos = build_road_graph(osm_nodes, osm_ways, from_wgs, risk_norm, gt)
    if not graph:
        print("[HATA] Graf boş.")
        return

    # ── OSM tesisler ──────────────────────────────────────────
    assembly_raw, safe_raw = fetch_osm_facilities(bbox_wgs84, from_wgs)

    # Güvenli toplanma noktaları (faz 1 hedefleri)
    assembly_nodes = filter_safe_nodes(assembly_raw, pos, risk_norm, gt,
                                       threshold=SAFE_RISK_THRESHOLD)
    # Güvenli acil tesisler (faz 2 hedefleri)
    safe_nodes     = filter_safe_nodes(safe_raw,     pos, risk_norm, gt,
                                       threshold=SAFE_RISK_THRESHOLD)

    # Fallback: toplanma noktası yoksa en yakın düşük riskli düğümler
    if not assembly_nodes:
        print("[Fallback] Toplanma noktası yok → düşük riskli yol düğümleri kullanılıyor.")
        candidates = [(nid, x, y)
                      for nid, (x, y) in pos.items()
                      if sample_risk(risk_norm, gt, x, y) < SAFE_RISK_THRESHOLD]
        # Uniform örnekleme
        step = max(1, len(candidates) // 50)
        for i in range(0, len(candidates), step):
            nid, x, y = candidates[i]
            assembly_nodes.append((nid, "low_risk_node", x, y))

    if not safe_nodes:
        print("[Fallback] Güvenli tesis yok → assembly_nodes hedef olarak kullanılıyor.")
        safe_nodes = assembly_nodes

    print(f"[Hedef] {len(assembly_nodes)} geçerli toplanma noktası, "
          f"{len(safe_nodes)} güvenli tesis")

    assembly_node_set = {n[0] for n in assembly_nodes}
    safe_node_set     = {n[0] for n in safe_nodes}

    # ── İki fazlı rotalama ────────────────────────────────────
    phase1_routes  = []   # halk → toplanma
    phase2_routes  = []   # devlet aracı → güvenli tesis
    assembly_used  = []   # QgsPoint layer için

    print(f"\n[Routing] {len(clusters)} küme işleniyor...")

    for idx, cluster in enumerate(clusters):
        cx, cy  = cluster_centroid_geo(cluster, gt)
        # Kaynak düğümü: küme merkezine en yakın yol düğümü
        src_nid, src_d = nearest_node(cx, cy, pos, max_dist=3000)
        if src_nid is None:
            print(f"  ✗ Küme {idx+1}: yakın yol düğümü yok (mesafe={src_d:.0f}m)")
            continue
        if src_nid not in graph:
            print(f"  ✗ Küme {idx+1}: kaynak düğüm grafta yok.")
            continue

        # ── Faz 1: Halk → Toplanma noktası ──────────────────
        path1, hit1 = dijkstra_graph(graph, [src_nid], assembly_node_set)
        if not path1:
            print(f"  ✗ Küme {idx+1}: toplanma noktasına rota bulunamadı")
            continue

        coords1    = path_nodes_to_coords(path1, pos)
        asm_label  = next((n[1] for n in assembly_nodes if n[0] == hit1), "toplanma")
        asm_x, asm_y = pos[hit1]

        phase1_routes.append({
            "id": idx, "phase": 1,
            "path_coords": coords1,
            "label": asm_label,
        })
        assembly_used.append((asm_x, asm_y, asm_label))

        print(f"  ✓ Küme {idx+1} Faz1 → {asm_label} ({len(coords1)} düğüm)")

        # ── Faz 2: Devlet aracı toplanma → güvenli tesis ────
        path2, hit2 = dijkstra_graph(graph, [hit1], safe_node_set)
        if path2:
            coords2   = path_nodes_to_coords(path2, pos)
            safe_lbl  = next((n[1] for n in safe_nodes if n[0] == hit2), "güvenli tesis")
            phase2_routes.append({
                "id": idx, "phase": 2,
                "path_coords": coords2,
                "label": safe_lbl,
            })
            print(f"           Faz2 → {safe_lbl} ({len(coords2)} düğüm)")
        else:
            print(f"           Faz2: güvenli tesise rota bulunamadı")

    # ── QGIS Katmanları ───────────────────────────────────────
    if not phase1_routes and not phase2_routes:
        print("\n[HATA] Hiçbir rota oluşturulamadı.")
        return

    proj = QgsProject.instance()

    if phase1_routes:
        lyr1 = make_route_layer(phase1_routes, "EPSG:25830",
                                "Faz 1 — Halk → Toplanma")
        style_route_layer(lyr1, phase=1)
        proj.addMapLayer(lyr1)

    if phase2_routes:
        lyr2 = make_route_layer(phase2_routes, "EPSG:25830",
                                "Faz 2 — Tahliye (Devlet Aracı)")
        style_route_layer(lyr2, phase=2)
        proj.addMapLayer(lyr2)

    if assembly_used:
        uniq_asm = list({(x, y, l) for x, y, l in assembly_used})
        lyr_asm  = make_point_layer(uniq_asm, "EPSG:25830", "Toplanma Noktaları")
        style_point_layer(lyr_asm, "0,200,100,255", size="5")
        proj.addMapLayer(lyr_asm)

    # Kaynak bölge merkezleri
    src_pts = []
    for i, cl in enumerate(clusters):
        cx, cy = cluster_centroid_geo(cl, gt)
        src_pts.append((cx, cy, f"risk_zone_{i+1}"))
    if src_pts:
        lyr_src = make_point_layer(src_pts, "EPSG:25830", "Risk Bölgeleri (Kaynak)")
        style_point_layer(lyr_src, "220,20,20,200", size="4")
        proj.addMapLayer(lyr_src)

    print("\n" + "=" * 65)
    print(f"✅ Tamamlandı!")
    print(f"   Faz 1 rotası : {len(phase1_routes)}  (turuncu kesik çizgi)")
    print(f"   Faz 2 rotası : {len(phase2_routes)}  (kırmızı düz çizgi)")
    print(f"   Toplanma noktaları: {len({(x,y) for x,y,_ in assembly_used})}")
    print(f"   Layer panelinde 4 yeni katman eklendi.")
    print("=" * 65)


# ── ÇALIŞTIR ──────────────────────────────────────────────────
run_evacuation_routing()
