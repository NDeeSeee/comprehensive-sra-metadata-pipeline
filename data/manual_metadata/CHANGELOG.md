  # CHANGELOG - POSEIDON Project

  ## [2024-12-19] - File Naming Standardization

  **CHANGES MADE:**
  - ✅ Renamed `oropharyngeal.xlsx` → `Oropharyngeal.xlsx` (case standardization)
  - ✅ Renamed `Adenocarcinoma of the Esophagus.xlsx` → `Esophagus.xlsx` (simplified naming)
  - ✅ Renamed `Gallbladder-SRA.xlsx` → `Gallbladder.xlsx` (removed SRA suffix)

  **REASONING:**
  - Standardized Excel file names to match existing directory structure
  - Ensures scripts can reliably find corresponding directories
  - Maintains consistency across the project
  - No data loss, only file renaming

  **IMPACT:**
  - All Excel files now have perfect 1:1 mapping with directory names
  - Scripts will work correctly with standardized naming
  - Improved maintainability and clarity

  ---

  **1. Excel Files vs Directory Names Comparison**

  | Excel File (SRAMetadataFiles/) | Directory Name | Status | Issue |
  |-----------------------------------|-------------------|------------|-----------|
  | Oropharyngeal.xlsx | Oropharyngeal/ | ✅ CONSISTENT | ✅ FIXED: Renamed from oropharyngeal.xlsx |
  | Esophagus.xlsx | Esophagus/ | ✅ CONSISTENT | ✅ FIXED: Renamed from "Adenocarcinoma of the Esophagus.xlsx" |
  | Gum+MouthOther.xlsx | Gum+MouthOther/ | ✅ CONSISTENT | - |
  | Gallbladder.xlsx | Gallbladder/ | ✅ CONSISTENT | ✅ FIXED: Renamed from "Gallbladder-SRA.xlsx" |
  | Larynx.xlsx | Larynx/ | ✅ CONSISTENT | - |
  | LungLargeCell.xlsx | LungLargeCell/ | ✅ CONSISTENT | - |
  | MouthFloor.xlsx | MouthFloor/ | ✅ CONSISTENT | - |
  | Pancreas.xlsx | Pancreas/ | ✅ CONSISTENT | - |
  | RenalPelvis.xlsx | RenalPelvis/ | ✅ CONSISTENT | - |
  | Salivary.xlsx | Salivary/ | ✅ CONSISTENT | - |
  | SmallCellLung.xlsx | SmallCellLung/ | ✅ CONSISTENT | - |
  | Tongue.xlsx | Tongue/ | ✅ CONSISTENT | - |
  | Vulva.xlsx | MISSING | ❌ MISSING | No directory exists |


  **2. Missing Cancer Types in Directories**

  | Excel File | Missing From | Impact |
  |---------------|------------------|------------|
  | Vulva.xlsx | All directories | No data processed yet |
  | Mesothelioma_controlled.csv | All directories | CSV format, needs conversion |
  | Pancreas_controlled.csv | All directories | CSV format, needs conversion |
  | Vagina_controls.csv | All directories | CSV format, needs conversion |
  | Vagina.csv | All directories | CSV format, needs conversion |

  **3. Directory Coverage Analysis**

  Tumors Directory (14 cancer types):
  • ✅ Complete coverage of most Excel files
  • ❌ Missing: Vulva, Mesothelioma (controlled), Vagina

  Controls Directory (9 cancer types):
  • ❌ Missing: Oropharyngeal, LungLargeCell, Mesothelioma, MouthFloor, RenalPelvis, Salivary, SmallCellLung, Tongue, Vulva

  Premalignant Directory (2 cancer types):
  • ❌ Very limited coverage: Only Esophagus, Mouth

  Bulk_CellTypes Directory (3 cancer types):
  • ❌ Limited coverage: Only Gallbladder, Mouth, Pancreas