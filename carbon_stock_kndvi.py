"""
Estimasi AGC dari AGB dengan kNDVI + Random Forest
====================================================
QGIS Processing Script

Letakkan file ini di:
  - Windows : %AppData%\QGIS\QGIS3\profiles\default\processing\scripts\
  - Linux/Mac: ~/.local/share/QGIS/QGIS3/profiles/default/processing/scripts\

Perbedaan dari versi NDVI:
  - Indeks vegetasi diganti kNDVI (kernel NDVI, RBF kernel)
  - Sigma dihitung pixel-wise: σ = 0.5 × (NIR + RED)
  - Fitur model RF: [RED, NIR, kNDVI]  (sama posisinya, beda indeksnya)

Referensi kNDVI:
  Camps-Valls et al. (2021). A unified vegetation index for quantifying
  the terrestrial biosphere. Science Advances, 7(9).
"""

from typing import Any, Optional
import numpy as np

from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingContext,
    QgsProcessingException,
    QgsProcessingFeedback,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterRasterDestination,
    QgsProcessingParameterString,
)
from qgis import processing

try:
    from osgeo import gdal
except ImportError:
    import gdal

try:
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import mean_squared_error, r2_score
except ImportError:
    raise QgsProcessingException("Library 'scikit-learn' belum terinstall di OSGeo4W.")


# ─────────────────────────────────────────────
# kNDVI Helper Functions
# ─────────────────────────────────────────────

def pixel_wise_sigma(nir: float, red: float) -> float:
    """
    Hitung sigma pixel-wise.
    σ = 0.5 × (NIR + RED)
    """
    return 0.5 * (nir + red)


def calculate_kndvi(nir: float, red: float) -> float:
    """
    Hitung kNDVI dengan RBF kernel.

    σ     = 0.5 × (NIR + RED)          ← pixel-wise
    knr   = exp(−(NIR − RED)² / (2σ²))
    kNDVI = (1 − knr) / (1 + knr)

    Args:
        nir : nilai piksel Near Infrared
        red : nilai piksel Red
    Returns:
        kndvi : float [−1, 1]
    """
    sigma = pixel_wise_sigma(nir, red)
    if sigma == 0:
        sigma = 1e-10
    knr   = np.exp(-((nir - red) ** 2) / (2 * sigma ** 2))
    return float((1 - knr) / (1 + knr))


def calculate_kndvi_array(nir_arr: np.ndarray, red_arr: np.ndarray) -> np.ndarray:
    """
    Versi vectorised untuk array raster satu baris.

    Args:
        nir_arr : array NIR (1D, satu baris piksel)
        red_arr : array RED (1D, satu baris piksel)
    Returns:
        kndvi   : array float [−1, 1], shape sama dengan input
    """
    sigma = 0.5 * (nir_arr + red_arr)
    sigma = np.where(sigma == 0, 1e-10, sigma)
    knr   = np.exp(-((nir_arr - red_arr) ** 2) / (2 * sigma ** 2))
    return (1 - knr) / (1 + knr)


# ─────────────────────────────────────────────
# Processing Algorithm
# ─────────────────────────────────────────────

class AGCCalculationAlgorithm(QgsProcessingAlgorithm):

    RED_BAND      = "RED_BAND"
    NIR_BAND      = "NIR_BAND"
    INPUT_VEKTOR  = "INPUT_VEKTOR"
    OUTPUT_RASTER = "OUTPUT_RASTER"

    # ── Metadata ──────────────────────────────
    def name(self) -> str:
        return "agc_from_agb_kndvi_model"

    def displayName(self) -> str:
        return "Estimasi AGC dari AGB (kNDVI + RF Model)"

    def group(self) -> str:
        return "Analisis Karbon"

    def groupId(self) -> str:
        return "analisiskarbon"

    def createInstance(self):
        return AGCCalculationAlgorithm()

    def shortHelpString(self) -> str:
        return (
            "<b>Estimasi AGC dari AGB menggunakan kNDVI + Random Forest</b><br><br>"
            "Script ini memprediksi peta Above Ground Carbon (AGC) dari data raster "
            "menggunakan model Random Forest yang dilatih dari titik sampel AGB lapangan.<br><br>"
            "<b>Fitur model RF:</b> RED, NIR, kNDVI<br>"
            "<b>kNDVI:</b> σ = 0.5×(NIR+RED) → kNDVI = (1−knr)/(1+knr)<br>"
            "<b>Konversi AGC:</b> AGC = AGB × 0.47 (IPCC 2006)<br><br>"
            "<b>Input:</b>"
            "<ul>"
            "  <li>Raster Band RED (misal Band 4 Sentinel-2)</li>"
            "  <li>Raster Band NIR (misal Band 8 Sentinel-2)</li>"
            "  <li>Layer titik sampel dengan field AGB (ton/ha)</li>"
            "</ul>"
            "<b>Output:</b> Raster prediksi AGC (Ton C/ha)<br><br>"
            "<b>Referensi:</b><br>"
            "Camps-Valls et al. (2021) — kNDVI, Science Advances<br>"
            "IPCC (2006) — National GHG Inventories"
        )

    # ── UI Parameters ─────────────────────────
    def initAlgorithm(self, config=None):

        # Band RED
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.RED_BAND,
                "Pilih Band RED (misal Band 4)"
            )
        )

        # Band NIR
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.NIR_BAND,
                "Pilih Band NIR (misal Band 8)"
            )
        )

        # Layer titik sampel AGB (ton/ha)
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_VEKTOR,
                "Input Titik Sampel AGB (ton/ha)",
                [QgsProcessing.SourceType.TypeVectorPoint]
            )
        )

        # Output raster
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT_RASTER,
                "Peta Prediksi AGC (Ton C/ha)"
            )
        )

    # ── Main Process ──────────────────────────
    def processAlgorithm(
        self,
        parameters: dict,
        context: QgsProcessingContext,
        feedback: QgsProcessingFeedback,
    ) -> dict:

        red_layer = self.parameterAsRasterLayer(parameters, self.RED_BAND, context)
        nir_layer = self.parameterAsRasterLayer(parameters, self.NIR_BAND, context)
        agb_field = "AGB"  # field AGB hardcode, satuan ton/ha

        # ── 1. Ekstraksi nilai piksel ke titik sampel ────────────────────
        feedback.pushInfo("=" * 50)
        feedback.pushInfo("LANGKAH 1: Ekstraksi nilai piksel ke titik sampel...")

        sampled_red = processing.run(
            "native:rastersampling",
            {
                "INPUT"        : parameters[self.INPUT_VEKTOR],
                "RASTERCOPY"   : parameters[self.RED_BAND],
                "COLUMN_PREFIX": "RED_",
                "OUTPUT"       : "memory:sampled_red",
            },
            context=context, feedback=feedback,
        )["OUTPUT"]

        sampled_all = processing.run(
            "native:rastersampling",
            {
                "INPUT"        : sampled_red,
                "RASTERCOPY"   : parameters[self.NIR_BAND],
                "COLUMN_PREFIX": "NIR_",
                "OUTPUT"       : "memory:sampled_all",
            },
            context=context, feedback=feedback,
        )["OUTPUT"]

        # ── 2. Persiapan data training ───────────────────────────────────
        feedback.pushInfo("LANGKAH 2: Persiapan data training...")
        X, y = [], []
        skipped = 0

        for feat in sampled_all.getFeatures():
            attrs    = feat.attributes()
            val_red  = attrs[-2]   # kolom RED_ (hasil sampling terakhir -2)
            val_nir  = attrs[-1]   # kolom NIR_ (hasil sampling terakhir -1)
            val_agb  = feat[agb_field]

            if None in [val_red, val_nir, val_agb]:
                skipped += 1
                continue

            red = float(val_red)
            nir = float(val_nir)

            # kNDVI — sigma pixel-wise
            kndvi = calculate_kndvi(nir, red)

            # Konversi AGB → AGC  (IPCC 2006: fraksi karbon = 0.47)
            agc = float(val_agb) * 0.47

            X.append([red, nir, kndvi])
            y.append(agc)

        X = np.array(X)
        y = np.array(y)

        feedback.pushInfo(f"  Sampel valid   : {len(y)}")
        feedback.pushInfo(f"  Sampel dilewati: {skipped}")

        if len(y) < 5:
            raise QgsProcessingException(
                f"Sampel terlalu sedikit ({len(y)} titik valid). "
                "Minimal 5 titik diperlukan untuk training RF."
            )

        # ── 3. Training Random Forest ────────────────────────────────────
        feedback.pushInfo("LANGKAH 3: Training model Random Forest...")
        feedback.pushInfo("  Fitur: [RED, NIR, kNDVI]")

        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.3, random_state=42
        )

        rf = RandomForestRegressor(n_estimators=100, random_state=42)
        rf.fit(X_train, y_train)

        y_pred_test = rf.predict(X_test)
        rmse = float(np.sqrt(mean_squared_error(y_test, y_pred_test)))
        r2   = float(r2_score(y_test, y_pred_test))

        feedback.pushInfo(f"  R²  : {r2:.4f}")
        feedback.pushInfo(f"  RMSE: {rmse:.4f} Ton C/ha")

        # ── 4. Prediksi raster per baris ─────────────────────────────────
        feedback.pushInfo("LANGKAH 4: Prediksi raster AGC...")

        ds_red = gdal.Open(red_layer.source())
        ds_nir = gdal.Open(nir_layer.source())

        rows         = ds_red.RasterYSize
        cols         = ds_red.RasterXSize
        geotransform = ds_red.GetGeoTransform()
        pixel_res_m  = abs(geotransform[1])                     # resolusi spasial (meter)
        pixel_area_ha = (pixel_res_m ** 2) / 10_000             # luas 1 piksel (ha)

        output_path = self.parameterAsOutputLayer(parameters, self.OUTPUT_RASTER, context)
        driver  = gdal.GetDriverByName("GTiff")
        out_ds  = driver.Create(output_path, cols, rows, 1, gdal.GDT_Float32)
        out_ds.SetGeoTransform(geotransform)
        out_ds.SetProjection(ds_red.GetProjection())

        total_agc_ton = 0.0
        total_area_ha = 0.0

        band_red = ds_red.GetRasterBand(1)
        band_nir = ds_nir.GetRasterBand(1)
        band_out = out_ds.GetRasterBand(1)

        for r in range(rows):
            if feedback.isCanceled():
                break

            red_row = band_red.ReadAsArray(0, r, cols, 1).ravel().astype(float)
            nir_row = band_nir.ReadAsArray(0, r, cols, 1).ravel().astype(float)

            # kNDVI vectorised per baris
            kndvi_row = calculate_kndvi_array(nir_row, red_row)
            kndvi_row = np.nan_to_num(kndvi_row)

            # Susun fitur [RED, NIR, kNDVI]
            row_X    = np.stack([red_row, nir_row, kndvi_row], axis=1)
            row_X    = np.nan_to_num(row_X)

            # Prediksi AGC (ton C/ha per piksel)
            row_pred = rf.predict(row_X)

            # Akumulasi statistik (piksel dengan AGC > 0)
            mask          = row_pred > 0
            total_agc_ton += float(np.sum(row_pred[mask]))
            total_area_ha += float(np.count_nonzero(mask)) * pixel_area_ha

            band_out.WriteArray(row_pred.reshape(1, cols), 0, r)
            feedback.setProgress(int((r / rows) * 100))

        # ── 5. Laporan hasil ─────────────────────────────────────────────
        avg_agc = total_agc_ton / total_area_ha if total_area_ha > 0 else 0.0

        feedback.pushInfo("\n" + "=" * 50)
        feedback.pushInfo("HASIL ANALISIS")
        feedback.pushInfo(f"  Indeks vegetasi    : kNDVI (RBF kernel, σ pixel-wise)")
        feedback.pushInfo(f"  Konversi AGC       : AGB × 0.47 (IPCC 2006)")
        feedback.pushInfo(f"  R²  model RF       : {r2:.4f}")
        feedback.pushInfo(f"  RMSE model RF      : {rmse:.4f} Ton C/ha")
        feedback.pushInfo("-" * 50)
        feedback.pushInfo(f"  Total luas (AGC>0) : {total_area_ha:.2f} Ha")
        feedback.pushInfo(f"  Total stok AGC     : {total_agc_ton:.2f} Ton C")
        feedback.pushInfo(f"  Rata-rata AGC/Ha   : {avg_agc:.2f} Ton C/Ha")
        feedback.pushInfo("=" * 50 + "\n")

        out_ds.FlushCache()
        out_ds = None

        return {self.OUTPUT_RASTER: output_path}
