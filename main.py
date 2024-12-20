import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from fpdf import FPDF
import pandas as pd
import numpy as np

# Function to read S-box from file
def read_sbox_from_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()
        sbox = []
        for line in lines:
            row = list(map(int, line.split()))
            sbox.append(row)
        return np.array(sbox)

# Function to calculate Hamming distance
def haming_distance(a, b):
    return np.sum(a != b)

# Convert S-box to binary form
def sbox_to_binary(Sbox):
    size = Sbox.shape
    bin_sbox = np.zeros((size[0], size[1], 8), dtype=int)
    for i in range(size[0]):
        for j in range(size[1]):
            bin_sbox[i, j] = [int(b) for b in format(Sbox[i, j], '08b')]
    return bin_sbox

# Walsh-Hadamard Transform
def walsh_hadamard_transform(vec):
    n = len(vec)
    H = np.array(vec)
    h = 1
    while h < n:
        for i in range(0, n, h * 2):
            for j in range(h):
                x = H[i + j]
                y = H[i + j + h]
                H[i + j] = x + y
                H[i + j + h] = x - y
        h *= 2
    return H

# Calculate Nonlinearity
def calculate_nonlinearity(Sbox):
    n = 8  # Input bit size
    m = 8  # Output bit size
    nl_values = []
    for output_bit in range(m):
        truth_table = []
        for row in range(Sbox.shape[0]):
            for col in range(Sbox.shape[1]):
                input_val = row * Sbox.shape[1] + col
                output_val = Sbox[row, col]
                output_bit_val = (output_val >> output_bit) & 1
                truth_table.append((-1) ** output_bit_val)
        # Transform Walsh-Hadamard
        wh = walsh_hadamard_transform(truth_table)
        max_wh = max(abs(wh))
        nl = (2 ** (n - 1)) - (max_wh / 2)
        nl_values.append(nl)
    return min(nl_values)

# Calculate SAC
def calculate_sac(Sbox_bin, n_bits):
    sac_total = 0
    total_comparisons = 0
    num_rows, num_cols, _ = Sbox_bin.shape

    for input_index in range(num_rows * num_cols):
        original_output = Sbox_bin[input_index // num_cols, input_index % num_cols]

        for bit_position in range(n_bits):
            flipped_input = np.copy(original_output)
            flipped_input[bit_position] ^= 1

            flipped_index = input_index ^ (1 << bit_position)
            flipped_index_row = flipped_index // num_cols
            flipped_index_col = flipped_index % num_cols
            flipped_output = Sbox_bin[flipped_index_row, flipped_index_col]

            bit_changes = haming_distance(original_output, flipped_output)
            sac_total += bit_changes / n_bits
            total_comparisons += 1

    sac_value = sac_total / total_comparisons
    return sac_value

# Function to calculate BIC-NL

def bic_nl(sbox):
    non_linearities = []
    for bit in range(8):
        outputs = np.array([((x >> bit) & 1) * 2 - 1 for x in sbox.flatten()])
        spectrum = walsh_hadamard_transform(outputs)
        nl = (256 - np.max(np.abs(spectrum))) // 2
        non_linearities.append(nl)
    return min(non_linearities)

# Function to calculate BIC-SAC
def bic_sac(sbox):
    n = len(sbox)
    bit_length = 8
    total_pairs = 0
    total_independence = 0

    for i in range(bit_length):
        for j in range(i + 1, bit_length):
            independence_sum = 0

            for x in range(n):
                for bit_to_flip in range(bit_length):
                    flipped_x = x ^ (1 << bit_to_flip)

                    y1 = sbox[x]
                    y2 = sbox[flipped_x]

                    b1_i = (y1 >> i) & 1
                    b1_j = (y1 >> j) & 1

                    b2_i = (y2 >> i) & 1
                    b2_j = (y2 >> j) & 1

                    independence_sum += ((b1_i ^ b2_i) ^ (b1_j ^ b2_j))

            pair_independence = independence_sum / (n * bit_length)
            total_independence += pair_independence
            total_pairs += 1

    bic_sac = total_independence / total_pairs
    return round(bic_sac, 5)

# Function LAP
def calculate_lap(sbox):
    size = len(sbox)
    lap_values = []

    for a in range(1, size):
        for b in range(1, size):
            count = 0
            for x in range(size):
                if bin(a & x).count('1') % 2 == bin(b & sbox[x]).count('1') % 2:
                    count += 1
            lap = abs(count - size / 2) / size
            lap_values.append(lap)

    return max(lap_values)

# Function DAP
def calculate_dap(sbox):
    size = len(sbox)
    dap_values = []

    for delta_in in range(1, size):
        for delta_out in range(size):
            count = 0
            for x in range(size):
                if sbox[x] ^ sbox[x ^ delta_in] == delta_out:
                    count += 1
            dap = count / size
            dap_values.append(dap)
    return max(dap_values)

# Function to select file and process S-box
def select_file_and_process():
    filepath = filedialog.askopenfilename(filetypes=[("Text files", "*.txt")])
    if not filepath:
        return

    global Sbox, bin_sbox, nonlinearity, sac, bic_nl_value, bic_sac_value, lap, dap

    # Read S-box from file
    Sbox = read_sbox_from_file(filepath)

    # Calculations
    bin_sbox = sbox_to_binary(Sbox)
    nonlinearity = calculate_nonlinearity(Sbox)
    sac = calculate_sac(bin_sbox, 8)
    bic_nl_value = bic_nl(Sbox.flatten())
    bic_sac_value = bic_sac(Sbox.flatten())
    lap = calculate_lap(Sbox.flatten())
    dap = calculate_dap(Sbox.flatten())

    calculate_and_display()

# Function for manual input
def manual_input():
    manual_input_window = tk.Toplevel(root)
    manual_input_window.title("Input Manual S-Box")
    manual_input_window.geometry("400x300")

    tk.Label(manual_input_window, text="Enter S-Box values (16x16 matrix, comma-separated):", wraplength=300).pack(pady=10)
    input_text = tk.Text(manual_input_window, height=10)
    input_text.pack(pady=10)

    def process_manual_input():
        manual_input_data = input_text.get("1.0", tk.END).strip()
        try:
            rows = manual_input_data.split("\n")
            sbox = np.array([list(map(int, row.replace(",", " ").split())) for row in rows])
            if sbox.shape == (16, 16):
                global Sbox, bin_sbox, nonlinearity, sac, bic_nl_value, bic_sac_value, lap, dap
                Sbox = sbox
                bin_sbox = sbox_to_binary(Sbox)
                nonlinearity = calculate_nonlinearity(Sbox)
                sac = calculate_sac(bin_sbox, 8)
                bic_nl_value = bic_nl(Sbox.flatten())
                bic_sac_value = bic_sac(Sbox.flatten())
                lap = calculate_lap(Sbox.flatten())
                dap = calculate_dap(Sbox.flatten())
                calculate_and_display()
                manual_input_window.destroy()
            else:
                messagebox.showerror("Error", "S-Box must be a 16x16 matrix.")
        except Exception as e:
            messagebox.showerror("Error", f"Error processing manual input: {e}")

    submit_button = ttk.Button(manual_input_window, text="Submit", command=process_manual_input)
    submit_button.pack(pady=10)

# Function to save results to PDF
def save_to_pdf():
    pdf = FPDF()
    pdf.add_page()
    
    # Title
    pdf.set_font("Arial", size=14)
    pdf.cell(200, 10, txt="Hasil Analisis S-Box", ln=True, align="C")
    pdf.ln(10)
    
    # Table
    pdf.set_font("Arial", size=12)
    
    # Define column widths
    col_width = pdf.w / 2.5  # Adjusted width
    th = 10  # Row height
    pdf.cell(col_width, th, "Metric", 1, 0, "C")
    pdf.cell(col_width, th, "Value", 1, 1, "C")

    results = {
        "Nonlinearity": nonlinearity,
        "SAC": sac,
        "BIC-NL": bic_nl_value,
        "BIC-SAC": bic_sac_value,
        "LAP": lap,
        "DAP": dap
    }

    for key, value in results.items():
        pdf.cell(col_width, th, key, 1, 0, "C")
        pdf.cell(col_width, th, f"{value:.5f}", 1, 1, "C")

    filepath = filedialog.asksaveasfilename(defaultextension=".pdf", filetypes=[("PDF files", "*.pdf")])
    if filepath:
        pdf.output(filepath)

# Function to save results to Excel
def save_to_excel():
    results = {
        "Metric": ["Nonlinearity", "SAC", "BIC-NL", "BIC-SAC", "LAP", "DAP"],
        "Value": [
            f"{nonlinearity:.5f}",
            f"{sac:.5f}",
            f"{bic_nl_value:.5f}",
            f"{bic_sac_value:.5f}",
            f"{lap:.5f}",
            f"{dap:.5f}"
        ]
    }

    df = pd.DataFrame(results)

    filepath = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel files", "*.xlsx")])
    if filepath:
        df.to_excel(filepath, index=False)

# Function to display calculation results
def calculate_and_display():
    results = {
        "Nonlinearity": nonlinearity,
        "SAC": sac,
        "BIC-NL": bic_nl_value,
        "BIC-SAC": bic_sac_value,
        "LAP": lap,
        "DAP": dap
    }

    # Clear previous table
    for row in tree.get_children():
        tree.delete(row)

    # Add data to table
    for key, value in results.items():
        tree.insert("", "end", values=(key, f"{value:.5f}"))

# Initialize application
root = tk.Tk()
root.title("📊 Analisis S-Box")
root.geometry("600x800")
root.configure(bg="#f0f4c3")  # Light green background

# Global style for ttk
style = ttk.Style()
style.theme_use("clam")  # Modern theme
style.configure("TLabel", background="#f0f4c3", font=("Arial", 12))
style.configure("TButton", font=("Arial", 12, "bold"), padding=10)
style.configure("Treeview", font=("Arial", 12), rowheight=30)
style.configure("Treeview.Heading", font=("Arial", 14, "bold"), background="#8bc34a", foreground="#ffffff")
style.map("TButton", background=[("active", "#ffcc80"), ("!active", "#ffb74d")])

# Header frame
header_frame = tk.Frame(root, bg="#8bc34a", pady=20)  # Dark green for header
header_frame.pack(fill="x")

# Main title label
title_label = ttk.Label(
    header_frame,
    text="📊 Hasil Analisis S-Box",
    font=("Arial", 22, "bold"),
    foreground="#ffffff",
    anchor="center"
)
title_label.pack()

# Table frame
table_frame = ttk.Frame(root, padding=20)
table_frame.pack(fill="both", expand=True)

# Treeview table
tree = ttk.Treeview(table_frame, columns=("Metric", "Value"), show="headings", height=12)
tree.heading("Metric", text="📈 Metric")
tree.heading("Value", text="📝 Value")
tree.column("Metric", anchor="center", width=250)
tree.column("Value", anchor="center", width=250)
tree.pack(fill="both", expand=True, padx=10, pady=10)

# Button frame
button_frame = tk.Frame(root, bg="#f0f4c3", pady=20)
button_frame.pack()

# File selection button
select_file_button = ttk.Button(
    button_frame,
    text="🗂️ Pilih File S-Box",
    command=select_file_and_process,
    style="TButton"
)
select_file_button.grid(row=0, column=0, padx=10, pady=10)

# Manual input button
manual_input_button = ttk.Button(
    button_frame,
    text="✏️ Input Manual S-Box",
    command=manual_input,
    style="TButton"
)
manual_input_button.grid(row=0, column=1, padx=10, pady=10)

# Save to PDF button
save_pdf_button = ttk.Button(
    button_frame,
    text="📄 Simpan ke PDF",
    command=save_to_pdf,
    style="TButton"
)
save_pdf_button.grid(row=0, column=2, padx=10, pady=10)

# Save to Excel button
save_excel_button = ttk.Button(
    button_frame,
    text="🧾 Simpan ke Excel",
    command=save_to_excel,
    style="TButton"
)
save_excel_button.grid(row=0, column=3, padx=10, pady=10)

# Footer frame
footer_frame = tk.Frame(root, bg="#ffcc80", pady=15)  # Light yellow for footer
footer_frame.pack(fill="x", side="bottom")

# Footer label
footer_label = ttk.Label(
    footer_frame,
    text="Oleh: Kelompok 11 (Reguler)\nChrisandito S. E. Bia - 4611422135, Maulana A. Ibrahim - 4611422140\nM. Haikal A. Juniartha - 4611422146, Naufal H. Muzakki - 4611422156",
    font=("Arial", 10),
    foreground="#000000",
    anchor="center"
)
footer_label.pack()

# Run the application
root.mainloop()
