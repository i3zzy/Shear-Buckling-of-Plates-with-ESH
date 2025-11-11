# Shear Buckling of Plates with Edge-Stiffened Holes

![Logo](Logo.png)
| Developed by **Abdulaziz Alabdulwahab**  
📧 Email: alabdulwahab01@gmail.com  
© 2025 — All Rights Reserved

---

## 📚 Description

This repository contains Python scripts developed for **Abaqus/CAE** to automate the modeling, meshing, and shear elastic buckling analysis of **flat plates with flanged (edge-stiffened) perforations**.  
It supports circular, square, and diamond holes with customizable stiffener geometries and includes an analytical–numerical comparison model to predict shear buckling coefficients and critical shear stresses. Analytical can be found here "##"

This work forms part of an ongoing research project at the **University of Sydney**, School of Civil Engineering.

---

## ⚙️ Features

- Fully parametric control over plate and hole geometry  
- Automatic part creation, stiffener extrusion, and meshing in Abaqus  
- Buckling step setup and analysis
- Compares analytical prediction model for shear buckling with FEM
- Supports **Circular**, **Square**, and **Diamond** hole shapes  

---

## 🧪 How to Use

1. Open **Abaqus/CAE**.  
2. Run the main file:  
File → Run Script → main_script.py
3. Enter geometric and meshing parameters in the popup window.  
4. After analysis, the script will automatically:
- Retrieve the first buckling eigenvalue  
- Compute the predicted and FEM-based `τ_cr` values  
- Print comparison results in the console

---

## 🧩 File Structure

| File | Description |
|------|--------------|
| `main_script.py` | Core Abaqus automation + shear buckling model |
| `LICENSE.txt` | License terms (educational / non-commercial) |
| `CITATION.cff` | Citation metadata for academic referencing |
| `.gitignore` | Ignores temp, ODB, and Abaqus cache files |
| `README.md` | This documentation file |

---

## 📝 License

This repository is released under the **Educational Use License (EUL)**.  
You are free to use, modify, and distribute the code **for research or educational purposes only**.  
Commercial, proprietary, or redistributive use without written permission from the author is strictly prohibited.  

See [`LICENSE.txt`](LICENSE.txt) for the full license text.

---

## 🧾 Citation

If you use this code or its methodology in your research, please cite it as:
Alabdulwahab, A. (2025). Shear Buckling of Plates with Edge-Stiffened Holes.
University of Sydney. GitHub Repository: [https://github.com/i3zzy/Shear-Buckling-of-Plates-with-ESH/](https://github.com/i3zzy/Shear-Buckling-of-Plates-with-ESH/)


---

## 💬 Contact

For questions, collaborations, or permissions:  
📧 **alabdulwahab01@gmail.com**

---

> *Developed with passion for advancing Cold-Formed Steel research.*  
> *Happy FE!*
