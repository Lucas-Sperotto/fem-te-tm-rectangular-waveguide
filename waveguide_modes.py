import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
import meshio
import subprocess
import pandas as pd
import os
from scipy.special import jn_zeros  # Importa função para obter raízes das funções de Bessel

# === Função para gerar a malha cretangular no GMSH ===
def generate_rectangle_mesh(a=0.02, b=0.01, filename='rect'):
    lc = min(a, b) / 20
    geo_content = f"""
    lc = {lc};
    Point(1) = {{0, 0, 0, lc}};
    Point(2) = {{{a}, 0, 0, lc}};
    Point(3) = {{{a}, {b}, 0, lc}};
    Point(4) = {{0, {b}, 0, lc}};
    Line(1) = {{1, 2}};
    Line(2) = {{2, 3}};
    Line(3) = {{3, 4}};
    Line(4) = {{4, 1}};
    Line Loop(5) = {{1, 2, 3, 4}};
    Plane Surface(6) = {{5}};
    Physical Curve(100) = {{1, 2, 3, 4}};
    Physical Surface(200) = {{6}};
    Mesh 2;
    Save "{filename}.msh";
    """
    with open(f"{filename}.geo", "w") as f:
        f.write(geo_content)
    subprocess.run(["gmsh", f"{filename}.geo", "-2", "-format", "msh2"], check=True)



# === Função principal para resolver os modos com GMSH e FEM ===
def solve_modes_with_gmsh(a=0.01, b=0.005, mode='TM', num_modes=6, filename='teste01'):
    a = float(np.squeeze(a))
    b = float(np.squeeze(b))
    # === Geração da malha com GMSH ===
    generate_rectangle_mesh(a, b, filename)
    mesh = meshio.read(f"{filename}.msh")

    # === Leitura dos nós e elementos ===
    points = mesh.points[:, :2]
    triangles = mesh.cells_dict["triangle"]
    Nn = len(points)  # Número de nós
    X0 = np.zeros(Nn)  # Vetor de solução

    # === Identificação dos nós da borda ===
    boundary_nodes = set()
    for cell_block, phys_ids in zip(mesh.cells, mesh.cell_data_dict["gmsh:physical"].values()):
        if cell_block.type == "line":
            for i, line in enumerate(cell_block.data):
                if phys_ids[i] == 100:
                    boundary_nodes.update(line)

    # === Modo TM: condição de Dirichlet (exclui nós da borda) ===
    if mode == 'TM':
        node_id = np.ones(Nn, dtype=int)
        for i in boundary_nodes:
            node_id[i] = 0
            X0[i] = 0

        # Indexação apenas dos nós internos (não-borda)
        index = np.zeros(Nn, dtype=int)
        counter = 0
        for i in range(Nn):
            if node_id[i] == 1:
                counter += 1
                index[i] = counter
        Nf = counter  # Graus de liberdade

        # Montagem das matrizes S e T com restrição de borda
        S = lil_matrix((Nf, Nf))
        T = lil_matrix((Nf, Nf))
        for tri in triangles:
            n = tri
            x = points[n, 0]
            y = points[n, 1]
            mat = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
            De = np.linalg.det(mat)
            Ae = abs(De / 2)
            b = np.array([(y[1] - y[2]), (y[2] - y[0]), (y[0] - y[1])])#b = np.array([(y[1] - y[2]) / De, (y[2] - y[0]) / De, (y[0] - y[1]) / De])
            c_ = np.array([(x[2] - x[1]), (x[0] - x[2]), (x[1] - x[0])])#c_ = np.array([(x[2] - x[1]) / De, (x[0] - x[2]) / De, (x[1] - x[0]) / De])
            for i in range(3):
                for j in range(3):
                    Se = (b[i] * b[j] + c_[i] * c_[j]) / (4 * Ae) #Se = (b[i] * b[j] + c_[i] * c_[j]) * Ae
                    Te = Ae / 6 if i == j else Ae / 12
                    if node_id[n[i]] != 0 and node_id[n[j]] != 0:
                        S[index[n[i]]-1, index[n[j]]-1] += Se
                        T[index[n[i]]-1, index[n[j]]-1] += Te

        # Autovalores apenas nos nós internos
        vals, vecs = eigsh(S, k=num_modes+4, M=T, sigma=0, which='LM')
        # === Exporta autovalores, raiz e kc * r ===
        autovalores_completos = np.real(vals)
        raiz_autovalores = np.sqrt(np.clip(autovalores_completos, 0, None))  # Garante não-negatividade
        kc_a = raiz_autovalores * a

        # Cria DataFrame com autovalores
        df_autovalores = pd.DataFrame({
            'Índice': list(range(1, len(autovalores_completos)+1)),
            'Autovalor (λ)': autovalores_completos,
            'Raiz (√λ = kc)': raiz_autovalores,
            'kc * r': kc_a
        })

        # Salva arquivo em CSV separado
        output_dir = f'out/results/{mode.lower()}_{filename}'
        os.makedirs(output_dir, exist_ok=True)
        autoval_path = os.path.join(output_dir, 'autovalores_detalhados.csv')
        df_autovalores.to_csv(autoval_path, index=False)

        kc = np.sqrt(np.real(vals[:num_modes]))
        fc = 3e8 * kc / (2 * np.pi)

        # Reconstrói a solução completa (incluindo zeros nos nós da borda)
        modos = []
        for q in range(num_modes):
            X0[:] = 0
            j = 0
            for i in range(Nn):
                if index[i] != 0:
                    X0[i] = vecs[j, q]
                    j += 1
            modos.append(X0.copy())

    # === Modo TE: condição de Neumann (todos os nós participam) ===
    elif mode == 'TE':
        # Todos os nós participam da solução (sem exclusão)
        S = lil_matrix((Nn, Nn))
        T = lil_matrix((Nn, Nn))
        for tri in triangles:
            n = tri
            x = points[n, 0]
            y = points[n, 1]
            mat = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
            De = np.linalg.det(mat)
            Ae = abs(De / 2)
            b = np.array([(y[1] - y[2]), (y[2] - y[0]), (y[0] - y[1])])#b = np.array([(y[1] - y[2]) / De, (y[2] - y[0]) / De, (y[0] - y[1]) / De])
            c_ = np.array([(x[2] - x[1]), (x[0] - x[2]), (x[1] - x[0])])#c_ = np.array([(x[2] - x[1]) / De, (x[0] - x[2]) / De, (x[1] - x[0]) / De])
            for i in range(3):
                for j in range(3):
                    Se = (b[i] * b[j] + c_[i] * c_[j]) / (4 * Ae) #Se = (b[i] * b[j] + c_[i] * c_[j]) * Ae
                    Te = Ae / 6 if i == j else Ae / 12
                    S[n[i], n[j]] += Se
                    T[n[i], n[j]] += Te

        # Autovalores no domínio completo (todos os nós)
        vals, vecs = eigsh(S, k=num_modes+4, M=T, sigma=0, which='LM')
        vals = np.real(vals)

        # === Exporta autovalores, raiz e kc * r ===
        autovalores_completos = np.real(vals)
        raiz_autovalores = np.sqrt(np.clip(autovalores_completos, 0, None))  # Garante não-negatividade
        kc_a = raiz_autovalores * a

        # Cria DataFrame com autovalores
        df_autovalores = pd.DataFrame({
            'Índice': list(range(1, len(autovalores_completos)+1)),
            'Autovalor (λ)': autovalores_completos,
            'Raiz (√λ = kc)': raiz_autovalores,
            'kc * r': kc_a
        })

        # Salva arquivo em CSV separado
        output_dir = f'out/results/{mode.lower()}_{filename}'
        os.makedirs(output_dir, exist_ok=True)
        autoval_path = os.path.join(output_dir, 'autovalores_detalhados.csv')
        df_autovalores.to_csv(autoval_path, index=False)

        # Elimina o primeiro autovalor (zero, associado ao modo trivial)
        kc = np.sqrt(vals[1:num_modes+1])
        fc = 3e8 * kc / (2 * np.pi)

        # Coleta os modos correspondentes
        modos = [vecs[:, i] for i in range(1, num_modes+1)]

    else:
        raise ValueError("Modo inválido. Use 'TM' ou 'TE'.")

    # === Diretório para salvar as imagens ===
    save_path = f"out/img/{mode.lower()}_{filename}"
    os.makedirs(save_path, exist_ok=True)

    # === Geração dos gráficos ===
    for q in range(num_modes):
        plt.figure(figsize=(6, 5))
        plt.tricontourf(points[:, 0], points[:, 1], triangles, modos[q], levels=100, cmap='jet')
        plt.colorbar(label='Amplitude')
         # Contorno do guia
        x_rect = [0.0, 0.01, 0.01, 0.0, 0.0]
        y_rect = [0.0, 0.0, 0.005, 0.005, 0.0]
        plt.plot(x_rect, y_rect, 'k--', linewidth=1)
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        plt.axis('equal')
        plt.title(f"Modo {q+1} - {mode}\nfc = {fc[q]/1e9:.3f} GHz | kc·a = {kc[q]*a:.3f}")
        plt.tight_layout()
        fig_path = os.path.join(save_path, f"modo_{q+1}_{mode}.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()

    # === Pasta para gráficos vetoriais ===
    save_path_quiver = f"out/img/quiver_{mode.lower()}_{filename}"
    os.makedirs(save_path_quiver, exist_ok=True)

    # === Calcula gradiente e gera quiver plot para cada modo ===
    from scipy.spatial import cKDTree

    tree = cKDTree(points)
    delta = a / 1000  # Passo pequeno para estimar derivadas

    # === Novo cálculo vetorial usando interpolação FEM (gradiente por triângulo) ===
    save_path_quiver = f"out/img/quiver_{mode.lower()}_{filename}"
    os.makedirs(save_path_quiver, exist_ok=True)

    from scipy.constants import mu_0, pi

    for q in range(num_modes):
        campo_z = modos[q]
        Ex = np.zeros(Nn)
        Ey = np.zeros(Nn)
        contagem = np.zeros(Nn)  # Contar quantas vezes cada nó participa (para média)

        for tri in triangles:
            n = tri
            x = points[n, 0]
            y = points[n, 1]
            mat = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
            De = np.linalg.det(mat)
            Ae = abs(De / 2)

            b = np.array([(y[1] - y[2]) / De, (y[2] - y[0]) / De, (y[0] - y[1]) / De])
            c_ = np.array([(x[2] - x[1]) / De, (x[0] - x[2]) / De, (x[1] - x[0]) / De])

            grad_phi_x = np.dot(campo_z[n], b)  # ∂φ/∂x
            grad_phi_y = np.dot(campo_z[n], c_)  # ∂φ/∂y

            if mode == 'TM':
                Ex_local = -grad_phi_y
                Ey_local = grad_phi_x
            elif mode == 'TE':
                kc_local = kc[q] if kc[q] != 0 else 1e-12
                omega = 2 * pi * fc[q]
                Z = omega * mu_0 / kc_local
                Ex_local = -Z * grad_phi_y
                Ey_local = -Z * grad_phi_x

            for i in n:
                Ex[i] += Ex_local
                Ey[i] += Ey_local
                contagem[i] += 1

        # Média dos valores nos nós
        Ex /= np.maximum(contagem, 1)
        Ey /= np.maximum(contagem, 1)

       
                # Cálculo da magnitude para normalização
        magnitude = np.sqrt(Ex**2 + Ey**2)
        nonzero = magnitude > 1e-14
        Ex_unit = np.zeros_like(Ex)
        Ey_unit = np.zeros_like(Ey)
        Ex_unit[nonzero] = Ex[nonzero] / magnitude[nonzero]
        Ey_unit[nonzero] = Ey[nonzero] / magnitude[nonzero]

        # Coloração pela magnitude original
        color_data = magnitude
        

                # === Cálculo da magnitude original dos vetores ===
        magnitude = np.sqrt(Ex**2 + Ey**2)
        max_mag = np.max(magnitude)

        # === Normalização de direção dos vetores ===
        Ex_unit = np.zeros_like(Ex)
        Ey_unit = np.zeros_like(Ey)
        nonzero = magnitude > 1e-14
        Ex_unit[nonzero] = Ex[nonzero] / magnitude[nonzero]
        Ey_unit[nonzero] = Ey[nonzero] / magnitude[nonzero]

        # === Ajuste de magnitude para visualização ===
        magnitude_norm = np.zeros_like(magnitude)
        magnitude_norm[nonzero] = magnitude[nonzero] / max_mag
        magnitude_scaled = 0.3 + 0.7 * magnitude_norm  # mínimo 30% do vetor máximo

        # === Vetores finais com tamanho ajustado ===
        Ex_final = Ex_unit * magnitude_scaled
        Ey_final = Ey_unit * magnitude_scaled

        # === Quiver Plot ===
        plt.figure(figsize=(6, 5))
        quiv = plt.quiver(
            points[:, 0], points[:, 1],
            Ex_final, Ey_final,
            magnitude,  # color map by original magnitude
            cmap='viridis',
            scale=30, width=0.002, edgecolor='k', linewidth=0.1
        )
        plt.colorbar(quiv, label='|Eₜ| (não normalizado)')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')

        # Contorno do guia
        x_rect = [0.0, 0.01, 0.01, 0.0, 0.0]
        y_rect = [0.0, 0.0, 0.005, 0.005, 0.0]
        plt.plot(x_rect, y_rect, 'k--', linewidth=1)

        plt.title(f"Campo Transversal - Modo {q+1} - {mode}\n"
                  f"fc = {fc[q]/1e9:.3f} GHz | kc·r = {kc[q]*a:.3f}")
        plt.axis('equal')
        plt.tight_layout()
        fig_path = os.path.join(save_path_quiver, f"quiver_modo_{q+1}_{mode}.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
    return fc, modos


# === Exporta os dados em CSV ===
def export_results_to_csv(fc, fcreal, error, kc, a, mode='TM', filename='teste01'):
    # Calcula kc * r (valor adimensional)
    kc_a = kc * a
    
    # Monta o dicionário com os dados
    data = {
        'Modo': [f'{i+1}' for i in range(len(fc))],
        'Frequência FEM (GHz)': fc / 1e9,
        'Frequência Teórica (GHz)': fcreal / 1e9,
        'kc * r': kc_a,
        'Erro Relativo (%)': error
    }

    # Cria DataFrame e salva em CSV
    df = pd.DataFrame(data)
    output_dir = f'out/results/{mode.lower()}_{filename}'
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, 'frequencias_modos.csv')
    df.to_csv(csv_path, index=False)
    return df


# === Bloco principal de execução ===
from scipy.special import jn_zeros

if __name__ == "__main__":
    # Parâmetros
    num_modos = 24
    b = 0.005  # altura do guia (m)
    a = 2 * b  # largura do guia (m)
    c = 3e8    # velocidade da luz

    # === Cálculo para modos TE ===
    print("Calculando modos TE...")
    fc_te, modos_te = solve_modes_with_gmsh(a=a, b=b, mode='TE', num_modes=num_modos, filename='te_retangular')

    # === Frequências de corte teóricas TE (m, n = 0, 1, 2, ...) ===
    modos_teoricos = []
    for m in range(0, 10):
        for n in range(1, 10):  # n ≠ 0 para TE
            if m == 0 and n == 0:
                continue
            fc = (c / 2) * np.sqrt((m / a)**2 + (n / b)**2)
            modos_teoricos.append((fc, m, n))
    modos_teoricos.sort()
    fcreal_TE = np.array([fc for fc, _, _ in modos_teoricos[:num_modos]])

    # Comparação
    kc_te_fem = 2 * np.pi * fc_te / c
    kc_te_teor = 2 * np.pi * fcreal_TE / c
    erro_TE = 100 * np.abs(fcreal_TE - fc_te) / fcreal_TE

    # Exportar resultados TE
    df_te = export_results_to_csv(fc_te, fcreal_TE, erro_TE, kc_te_fem, 1.0, mode='TE', filename='te_24modos')  # radius não se aplica, use 1.0 como dummy
    print("Resultados dos modos TE:")
    print(df_te)

    # === Cálculo para modos TM ===
    print("\nCalculando modos TM...")
    fc_tm, modos_tm = solve_modes_with_gmsh(a=a, b=b, mode='TM', num_modes=num_modos, filename='tm_24modos')

    # === Frequências de corte teóricas TM (m, n = 1, 2, ...) ===
    modos_teoricos = []
    for m in range(1, 10):
        for n in range(1, 10):
            fc = (c / 2) * np.sqrt((m / a)**2 + (n / b)**2)
            modos_teoricos.append((fc, m, n))
    modos_teoricos.sort()
    fcreal_TM = np.array([fc for fc, _, _ in modos_teoricos[:num_modos]])

    kc_tm_fem = 2 * np.pi * fc_tm / c
    kc_tm_teor = 2 * np.pi * fcreal_TM / c
    erro_TM = 100 * np.abs(fcreal_TM - fc_tm) / fcreal_TM

    # Exportar resultados TM
    df_tm = export_results_to_csv(fc_tm, fcreal_TM, erro_TM, kc_tm_fem, 1.0, mode='TM', filename='tm_24modos')
    print("Resultados dos modos TM:")
    print(df_tm)
