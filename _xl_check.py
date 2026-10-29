try:
    import openpyxl; print("openpyxl", openpyxl.__version__)
except Exception as e: print("no openpyxl:", e)
try:
    import pandas; print("pandas", pandas.__version__)
except Exception as e: print("no pandas:", e)
