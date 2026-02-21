# routeevacuationfromfloodriskmap
Flood Evacuation Routing — Two-Phase Road Network Optimization in QGIS

This project develops a GIS-based evacuation routing framework for flood disaster scenarios in Valencia, Spain. The model integrates raster-based population exposure, flood risk analysis, and OpenStreetMap road networks to compute optimal evacuation paths using a risk-weighted Dijkstra algorithm.

📌 Project Overview

Flood disasters require structured evacuation planning that accounts for:

Population density

Spatial flood risk

Road network constraints

Safe assembly and emergency facilities

This project implements a two-phase evacuation strategy:

Phase 1 — Civilian Evacuation

High-risk population clusters → Road network → Safe assembly points

Phase 2 — Government Evacuation Transport

Assembly points → Road network → Secure emergency facilities (hospital, police, fire station, shelter)

Flood risk is directly incorporated into routing cost functions to avoid high-risk road segments.

🗺 Study Area

Valencia, Spain
Coordinate Reference System: EPSG:25830 (ETRS89 / UTM zone 30N)

⚙️ Methodology
1️⃣ Risk–Population Overlap Analysis

HRSL population raster

Flood risk raster

Normalization & percentile thresholding

High-risk clusters detected via BFS region grouping

2️⃣ Road Network Extraction

OpenStreetMap data retrieved via Overpass API

Motorways, primary, secondary, residential, pedestrian paths included

One-way restrictions respected

3️⃣ Risk-Aware Graph Construction

Edge cost formula:

Cost = Distance × RoadSpeedFactor × (1 + FloodWeight × AverageFloodRisk)

This ensures routes avoid high flood-risk segments.

4️⃣ Facility Detection

OSM facilities classified as:

Assembly Points

emergency=assembly_point

parks

schools

stadiums

parking areas

Safe Facilities

hospitals

fire stations

police stations

clinics

shelters

Flood risk threshold filtering ensures safe destinations.

5️⃣ Routing Algorithm

Weighted Dijkstra shortest path

Multi-source to multi-destination routing

Two-phase sequential optimization

🧠 Algorithm Workflow

Load raster layers

Identify high-risk population clusters

Extract OSM road network

Build weighted graph

Detect safe assembly & emergency facilities

Run Phase 1 routing

Run Phase 2 routing

Generate styled QGIS output layers

📦 Requirements

QGIS 3.x

Python (bundled with QGIS)

NumPy

GDAL

requests

Internet connection (for Overpass API)

📂 Required Input Layers

Must be loaded in QGIS before running:

valenciahrsl → Population raster

valenciariskflood → Flood risk raster

Both in EPSG:25830
