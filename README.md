# Estimasi AGC dari AGB — kNDVI + Random Forest (QGIS Processing Script)

Script QGIS untuk memprediksi peta **Above Ground Carbon (AGC)** dari titik sampel biomassa lapangan (**AGB**) menggunakan indeks vegetasi **kNDVI** (kernel NDVI) dan model **Random Forest**.

---

## Deskripsi

Script ini merupakan QGIS Processing Algorithm yang mengotomatiskan alur estimasi stok karbon hutan dari data penginderaan jauh. Berbeda dari pendekatan konvensional yang menggunakan NDVI, script ini menggunakan **kNDVI berbasis RBF kernel** yang terbukti lebih stabil dan tidak jenuh pada vegetasi rapat (Camps-Valls et al., 2021).

### Alur Kerja

```
Band RED + Band NIR
        │
        ▼
Ekstraksi nilai piksel ke titik sampel (raster sampling)
        │
        ▼
Hitung kNDVI per titik sampel
σ = 0.5 × (NIR + RED)
knr = exp(−(NIR−RED)² / (2σ²))
kNDVI = (1 − knr) / (1 + knr)
        │
        ▼
Konversi AGB → AGC
AGC = AGB × 0.47  (IPCC 2006)
        │
        ▼
Training Random Forest
Fitur: [RED, NIR, kNDVI]
Target: AGC (ton C/ha)
        │
        ▼
Prediksi raster AGC seluruh citra
        │
        ▼
Output: Raster AGC (Ton C/ha) + Laporan statistik
```

---

## Rumus

### kNDVI (kernel NDVI — RBF kernel)

$$\sigma = 0.5 \cdot (NIR + RED)$$

$$k_{nr} = \exp\left(-\frac{(NIR - RED)^2}{2\sigma^2}\right)$$

$$kNDVI = \frac{1 - k_{nr}}{1 + k_{nr}}$$

### Konversi AGB ke AGC (IPCC 2006)

$$AGC = AGB \times 0.47$$

$$C_{BGB} = C_{AGB} \times 0.26$$

$$C_{total} = C_{AGB} + C_{BGB}$$

---

## Persyaratan

| Komponen | Versi |
|---|---|
| QGIS | ≥ 3.16 |
| Python | ≥ 3.8 (bawaan QGIS) |
| numpy | bawaan OSGeo4W |
| gdal / osgeo | bawaan OSGeo4W |
| scikit-learn | perlu diinstall manual |

### Install scikit-learn di OSGeo4W

Buka **OSGeo4W Shell** lalu jalankan:

```bash
pip install scikit-learn
```

---

## Instalasi Script

1. Unduh file `carbon_stock_kndvi.py`
2. Salin ke folder scripts QGIS:
   - **Windows**: `%AppData%\QGIS\QGIS3\profiles\default\processing\scripts\`
   - **Linux/Mac**: `~/.local/share/QGIS/QGIS3/profiles/default/processing/scripts/`
3. Buka QGIS → **Processing Toolbox** → klik ikon refresh ↺
4. Cari di: **Scripts → Analisis Karbon → Estimasi AGC dari AGB (kNDVI + RF Model)**

---

## Parameter Input

| Parameter | Tipe | Keterangan |
|---|---|---|
| Pilih Band RED | Raster Layer | Band merah (misal Band 4 Sentinel-2) |
| Pilih Band NIR | Raster Layer | Band inframerah dekat (misal Band 8 Sentinel-2) |
| Input Titik Sampel AGB (ton/ha) | Vector Point | Layer titik dengan field **AGB** (wajib bernama `AGB`, satuan ton/ha) |

> **Penting:** Kolom AGB di data vektor harus bernama `AGB` (huruf kapital semua).

---

## Output

| Output | Tipe | Keterangan |
|---|---|---|
| Peta Prediksi AGC | Raster GeoTIFF | Nilai AGC per piksel dalam satuan **Ton C/ha** |

### Laporan di Log QGIS

Setelah proses selesai, log akan menampilkan:

```
==================================================
HASIL ANALISIS
  Indeks vegetasi    : kNDVI (RBF kernel, σ pixel-wise)
  Konversi AGC       : AGB × 0.47 (IPCC 2006)
  R²  model RF       : 0.xxxx
  RMSE model RF      : x.xxxx Ton C/ha
--------------------------------------------------
  Total luas (AGC>0) : xxx.xx Ha
  Total stok AGC     : xxx.xx Ton C
  Rata-rata AGC/Ha   : xx.xx Ton C/Ha
==================================================
```

---

## Struktur Data Vektor

Contoh format atribut layer titik sampel yang diperlukan:

| FID | AGB |
|---|---|
| 1 | 120.5 |
| 2 | 85.3 |
| 3 | 210.0 |

Koordinat titik harus berada di dalam extent raster NIR dan RED.

---

## Referensi

- Camps-Valls, G., Campos-Taberner, M., Moreno-Martínez, Á., Walther, S., Duveiller, G., Cescatti, A., ... & Running, S. W. (2021). A unified vegetation index for quantifying the terrestrial biosphere. *Science Advances*, 7(9), eabc7447.
- IPCC. (2006). *2006 IPCC Guidelines for National Greenhouse Gas Inventories, Volume 4: Agriculture, Forestry and Other Land Use*. Chapter 4, Table 4.3 & 4.4.

---

## Lisensi

MIT License — bebas digunakan dan dimodifikasi dengan menyertakan atribusi.
