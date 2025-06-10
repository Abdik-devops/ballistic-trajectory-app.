import math
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from io import BytesIO
from openpyxl import Workbook

# Физические константы
R = 287.05287  # Газовая постоянная, Дж/(кг·К)
sa = 110.4     # Постоянная Сатерленда, К
betas = 1.458e-6  # Коэффициент вязкости, кг/(м·с·K^0.5)
gamma_air = 1.4  # Показатель адиабаты

def atm(h):
    """Расчет параметров атмосферы на заданной высоте"""
    try:
        if h < -2000 or h > 85000:
            raise ValueError(f"Высота {h} м вне допустимого диапазона [-2000, 85000] м")
        hg = Rz * h / (Rz + h)
        H = [-2000, 0, 11000, 20000, 32000, 47000, 51000, 71000, 85000]
        T = [301.15, 288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.65]
        P = [127773.7, 101325.0, 22632.1, 5474.877, 868.0158, 110.9058, 66.939, 3.9564, 0.3734]
        BETA = [-0.0065, -0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002, 0.0]
        for j in range(len(H)-1):
            if H[j] <= hg < H[j+1]:
                t = T[j] + BETA[j] * (hg - H[j])
                if BETA[j] != 0:
                    p = P[j] * (1 + BETA[j] * (hg - H[j]) / T[j])**(-gs / (R * BETA[j]))
                else:
                    p = P[j] * math.exp(-gs * (hg - H[j]) / (R * T[j]))
                break
        ro = p / (R * t)
        a = math.sqrt(gamma_air * R * t)
        mu = betas * t**1.5 / (t + sa)
        nu = mu / ro
        return t, p, ro, mu, nu, a
    except Exception as e:
        st.error(f"Ошибка в atm(): {e}")
        return 216.65, 22632.1, 0.3639, 1.789e-5, 4.917e-5, 295.1

def adh(M, q, omegax, omegay, omegaz, alpha, beta, deltav, deltan, deltae):
    """Расчет аэродинамических сил и моментов"""
    try:
        M = max(min(M, 20), 0.1)
        alpha = max(min(alpha, math.radians(89)), math.radians(-89))
        alpha_deg = math.degrees(alpha)
        Ds = 1.86 * (11.554 / math.exp(min(M, 10)) - 2.5191e-3 * M**2 - 5.024 / M + 52.836e-3 * M + 4.112)
        cyalpha = Ds**0.5 if Ds >= 0 else 1.86 * 1.039
        p1 = 1 / (243.84e-3 / math.exp(min(alpha_deg, 50)) + 74.309e-3)
        log_arg = max(1.9773 * alpha_deg**2 - 25.587 * alpha_deg + 83.354, 1e-10)
        p2 = math.log(log_arg)
        p3 = 18.985 * alpha_deg**2 - 375.76 * alpha_deg + 1471
        p4 = -51.164e-3 * alpha_deg**2 + 805.52e-3 * alpha_deg + 1.8929
        cydeltav = (-p1 * 1e-6 * M**2 + p2 * 1e-12 * math.exp(min(M, 10)) - p3 * 1e-6 * M - p4 * 1e-3) * 2
        czbeta = -cyalpha
        czdeltan = -cydeltav
        cx = 1 / (73.211 / math.exp(min(M, 10)) - 47.483 / M + 16.878)
        cy = cyalpha * alpha + cydeltav * deltav
        cz = czbeta * beta + czdeltan * deltan
        X = cx * q * s
        Y = cy * q * s
        Z = cz * q * s
        mxomegax = -0.005 * 0.6786
        mzomegaz = 1.89 * (146.79e-6 * M**2 - 158.98e-3 / M - 7.639e-3 * M - 68.195e-3)
        mzalpha = -766.79e-3 / math.exp(min(M, 10)) + 438.74e-3 / M + 5.8822e-3 * M - 158.34e-3
        k1 = math.exp(-19.488e-3 * alpha_deg**2 - 378.62e-3 * alpha_deg + 6.7518)
        k2 = math.exp(-21.234e-3 * alpha_deg**2 - 635.84e-6 * math.exp(min(alpha_deg, 50)) - 98.296e-3 * alpha_deg + 2.5938)
        mzdeltav = 1.89 * math.sqrt(max(k1 * 1e-9 * M**2 + k2 * 1e-6, 0))
        myomegay = mzomegaz
        mydeltan = mzdeltav
        mybeta = mzalpha
        Mstab = 0 if abs(omegax) >= 0.001 else 100
        v_safe = max(v, 1e-10)
        mx = mxomegax * omegax * l / v_safe
        my = myomegay * omegay * l / v_safe + mybeta * beta + mydeltan * deltan
        mz = mzomegaz * omegaz * l / v_safe + mzalpha * alpha + mzdeltav * deltav
        Mx = mx * q * s * l + Mstab * deltae
        My = my * q * s * l
        Mz = mz * q * s * l
        return X, Y, Z, Mx, My, Mz
    except Exception as e:
        st.error(f"Ошибка в adh(): {e}")
        return 0, 0, 0, 0, 0, 0

def matr(rorg, larg, murg, nurg):
    """Построение матрицы поворота"""
    m00 = rorg**2 + larg**2 - murg**2 - nurg**2
    m01 = 2 * (rorg * nurg + larg * murg)
    m02 = 2 * (-rorg * murg + larg * nurg)
    m10 = 2 * (-rorg * nurg + larg * murg)
    m11 = rorg**2 + murg**2 - nurg**2 - larg**2
    m12 = 2 * (rorg * larg + nurg * murg)
    m20 = 2 * (rorg * murg + nurg * larg)
    m21 = 2 * (-rorg * larg + nurg * murg)
    m22 = rorg**2 + nurg**2 - larg**2 - murg**2
    return np.array([[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]])

def eiler(d, b, c, e, dd, db, dc, de, dt):
    """Интегрирование методом Эйлера"""
    d += dd * dt
    b += db * dt
    c += dc * dt
    e += de * dt
    return d, b, c, e

def norm(ro, la, mu, nu):
    """Нормализация кватернионов"""
    u = max((ro**2 + la**2 + mu**2 + nu**2)**0.5, 1e-10)
    return ro / u, la / u, mu / u, nu / u

# Streamlit интерфейс
st.title("Симуляция баллистической траектории")
st.write("Введите параметры и нажмите 'Запустить симуляцию' для получения результатов.")

# Входные параметры
st.header("Входные параметры")
col1, col2 = st.columns(2)
with col1:
    N = st.number_input("N", value=9, min_value=0, max_value=100)
    theta_deg = st.number_input("Угол θ (градусы)", value=-41.0, min_value=-89.0, max_value=89.0)
    m = st.number_input("Масса (кг)", value=1242.0, min_value=1.0)
    v = st.number_input("Начальная скорость (м/с)", value=float(4850 - 5 * N), min_value=0.0)
    dt = st.number_input("Шаг времени (с)", value=0.001, min_value=1e-6, format="%.6f")
with col2:
    dm = st.number_input("Диаметр (м)", value=1.41, min_value=0.01)
    y0 = st.number_input("Начальная высота (м)", value=float(41000 - 52 * N), min_value=0.0)
    K1 = st.number_input("K1", value=-7.0)
    K2 = st.number_input("K2", value=-7.0)
    max_iterations = st.number_input("Макс. итераций", value=1000000, min_value=1000)

# Кнопка для запуска симуляции
if st.button("Запустить симуляцию"):
    with st.spinner("Выполняется симуляция..."):
        # Инициализация параметров
        theta = math.radians(theta_deg)
        Ix = 180
        Iy = Iz = 700
        l = 3.8
        deltav = deltan = deltae = 0
        K3 = -7
        K4 = -7
        wxg = 5
        wyg = wzg = 0
        psi = gamma = alpha = beta = omegax = omegay = omegaz = 0
        x = z = 0
        y = y0
        Rz = 6356767
        gs = 9.80665
        s = math.pi * dm**2 / 4
        Gx = Gz = 0
        Gy = -m * gs
        thetas = theta - alpha
        vy = v * math.sin(thetas)
        vx = v * math.cos(thetas)
        vz = 0
        vremya = 0

        # Массивы для хранения данных
        t_arr = [vremya]
        x_arr = [x]
        y_arr = [y]
        z_arr = [z]
        theta_arr = [math.degrees(theta)]
        psi_arr = [math.degrees(psi)]
        gamma_arr = [math.degrees(gamma)]
        alpha_arr = [math.degrees(alpha)]
        beta_arr = [math.degrees(beta)]
        vx_arr = [vx]
        vy_arr = [vy]
        vz_arr = [vz]
        deltav_arr = [math.degrees(deltav)]
        deltan_arr = [math.degrees(deltan)]
        deltae_arr = [math.degrees(deltae)]

        # Кватернионы
        rorg = math.cos(psi/2) * math.cos(theta/2) * math.cos(gamma/2) - math.sin(psi/2) * math.sin(theta/2) * math.sin(gamma/2)
        larg = math.sin(psi/2) * math.sin(theta/2) * math.cos(gamma/2) + math.cos(psi/2) * math.cos(theta/2) * math.sin(gamma/2)
        murg = math.sin(psi/2) * math.cos(theta/2) * math.cos(gamma/2) + math.cos(psi/2) * math.sin(theta/2) * math.sin(gamma/2)
        nurg = math.cos(psi/2) * math.sin(theta/2) * math.cos(gamma/2) - math.sin(psi/2) * math.cos(theta/2) * math.sin(gamma/2)

        # Основной цикл
        num = 0
        while y >= 0 and num < max_iterations:
            try:
                if y > 1e6 or v > 1e5:
                    st.error(f"Прерывание: y={y:.2f} м или v={v:.2f} м/с на итерации {num}")
                    break
                A = matr(rorg, larg, murg, nurg)
                det_A = np.linalg.det(A)
                if abs(det_A) < 1e-10:
                    st.error(f"Матрица A плохо обусловлена на итерации {num}, det(A)={det_A}")
                    break
                wg = np.array([[wxg], [wyg], [wzg]])
                wsv = np.dot(A, wg)
                wx = wsv[0, 0]
                wy = wsv[1, 0]
                wz = wsv[2, 0]
                v = ((vx + wx)**2 + (vy + wy)**2 + (vz + wz)**2)**0.5
                t, p, ro, mu, nu, a = atm(y)
                M = v / max(a, 1e-10)
                q = ro * v**2 / 2
                X, Y, Z, Mx, My, Mz = adh(M, q, omegax, omegay, omegaz, alpha, beta, deltav, deltan, deltae)
                Fsv = np.array([[-X], [Y], [Z]])
                F = np.dot(np.linalg.inv(A), Fsv)
                Fx = F[0, 0]
                Fy = F[1, 0]
                Fz = F[2, 0]
                drorgdt = -(omegax * larg + omegay * murg + omegaz * nurg) / 2
                dlargdt = (omegax * rorg - omegay * nurg + omegaz * murg) / 2
                dmurgdt = (omegax * nurg + omegay * rorg - omegaz * larg) / 2
                dnurgdt = (-omegax * murg + omegay * larg + omegaz * rorg) / 2
                dvxdt = (Fx + Gx) / m
                dvydt = (Fy + Gy) / m
                dvzdt = (Fz + Gz) / m
                domegaxdt = Mx / Ix - (Iz - Iy) * omegay * omegaz / Ix
                domegaydt = My / Iy - (Ix - Iz) * omegax * omegaz / Iy
                domegazdt = Mz / Iz - (Iy - Ix) * omegax * omegay / Iz
                rorg, larg, murg, nurg = eiler(rorg, larg, murg, nurg, drorgdt, dlargdt, dmurgdt, dnurgdt, dt)
                rorg, larg, murg, nurg = norm(rorg, larg, murg, nurg)
                vx, vy, vz = eiler(vx, vy, vz, 0, dvxdt, dvydt, dvzdt, 0, dt)[:3]
                vy = max(min(vy, 1e4), -1e4)
                omegax, omegay, omegaz = eiler(omegax, omegay, omegaz, 0, domegaxdt, domegaydt, domegazdt, 0, dt)[:3]
                x, y, z = eiler(x, y, z, 0, vx, vy, vz, 0, dt)[:3]
                denom_psi = rorg**2 + larg**2 - murg**2 - nurg**2
                denom_gamma = rorg**2 + murg**2 - nurg**2 - larg**2
                theta = math.asin(max(min(2 * (rorg * nurg + larg * murg), 1), -1))
                gamma = math.atan(2 * (rorg * larg - nurg * murg) / (denom_gamma if abs(denom_gamma) > 1e-10 else 1e-10))
                psi = math.atan(2 * (rorg * murg - larg * nurg) / (denom_psi if abs(denom_psi) > 1e-10 else 1e-10))
                dthetadt = omegay * math.sin(gamma) + omegaz * math.cos(gamma)
                dpsidt = (omegay * math.cos(gamma) - omegaz * math.sin(gamma)) / max(abs(math.cos(theta)), 1e-10)
                dgammadt = omegax - math.tan(theta) * (omegay * math.cos(gamma) - omegaz * math.sin(gamma))
                V = np.array([[vx], [vy], [vz]])
                Vsv = np.dot(A, V)
                vxsv = Vsv[0, 0]
                vysv = Vsv[1, 0]
                vzsv = Vsv[2, 0]
                alpha = -math.atan2(vysv + wy, max(vxsv + wx, 1e-10))
                v_total = math.sqrt(max((vxsv + wx)**2 + (vysv + wy)**2 + (vzsv + wz)**2, 1e-20))
                beta = math.asin((vzsv + wz) / v_total)
                deltav = -K1 * dthetadt
                deltan = -K2 * dpsidt
                deltae = -K3 * gamma - K4 * dgammadt
                vremya += dt

                if num % int(1 / dt) == 0:
                    t_arr.append(vremya)
                    x_arr.append(x)
                    y_arr.append(y)
                    z_arr.append(z)
                    theta_arr.append(math.degrees(theta))
                    psi_arr.append(math.degrees(psi))
                    gamma_arr.append(math.degrees(gamma))
                    alpha_arr.append(math.degrees(alpha))
                    beta_arr.append(math.degrees(beta))
                    vx_arr.append(vx)
                    vy_arr.append(vy)
                    vz_arr.append(vz)
                    deltav_arr.append(math.degrees(deltav))
                    deltan_arr.append(math.degrees(deltan))
                    deltae_arr.append(math.degrees(deltae))
                num += 1
            except Exception as e:
                st.error(f"Ошибка на итерации {num}: {e}")
                break

        # Создание DataFrame
        data = {
            'Время, с': t_arr,
            'Х, м': x_arr,
            'Y, м': y_arr,
            'Z, м': z_arr,
            'vx, м/с': vx_arr,
            'vy, м/с': vy_arr,
            'vz, м/с': vz_arr,
            'Тета, град': theta_arr,
            'Пси, град': psi_arr,
            'Гамма, град': gamma_arr,
            'Aльфа, град': alpha_arr,
            'Бета, град': beta_arr,
            'Дельта по тангажу, град': deltav_arr,
            'Дельта по рысканью, град': deltan_arr,
            'Дельта по крену, град': deltae_arr
        }
        df = pd.DataFrame(data)

        # Вывод результатов
        st.header("Результаты симуляции")
        st.subheader("Таблица данных")
        st.dataframe(df)

        # Графики
        st.subheader("Графики")
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        ax1.plot(df['Х, м'], df['Y, м'], label='Траектория (x-y)')
        ax1.set_xlabel('X, м')
        ax1.set_ylabel('Y, м')
        ax1.legend()
        ax1.grid(True)
        ax2.plot(df['Время, с'], df['Тета, град'], label='Тета')
        ax2.plot(df['Время, с'], df['Пси, град'], label='Пси')
        ax2.plot(df['Время, с'], df['Гамма, град'], label='Гамма')
        ax2.set_xlabel('Время, с')
        ax2.set_ylabel('Углы, град')
        ax2.legend()
        ax2.grid(True)
        plt.tight_layout()
        st.pyplot(fig)

        # Сохранение в Excel
        st.subheader("Скачать результаты")
        output = BytesIO()
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Параметры', index=False)
        excel_data = output.getvalue()
        st.download_button(
            label="Скачать Excel-файл",
            data=excel_data,
            file_name="UTS.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )

        st.success("Симуляция завершена!")