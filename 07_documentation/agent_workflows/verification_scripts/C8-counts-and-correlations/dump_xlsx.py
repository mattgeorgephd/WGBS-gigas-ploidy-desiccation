import sys, openpyxl, warnings
warnings.filterwarnings("ignore")
path = sys.argv[1]
maxrows = int(sys.argv[2]) if len(sys.argv) > 2 else 60
for data_only in (True,):
    wb = openpyxl.load_workbook(path, data_only=data_only, read_only=False)
    for ws in wb.worksheets:
        print(f"=== sheet '{ws.title}' state={ws.sheet_state} dims={ws.dimensions} max_row={ws.max_row} max_col={ws.max_column}")
        n = 0
        for row in ws.iter_rows():
            vals = [(c.coordinate, c.value) for c in row if c.value is not None]
            if vals:
                print("  " + " | ".join(f"{k}={v!r}" for k, v in vals)[:600])
                n += 1
                if n >= maxrows:
                    print("  ..."); break
