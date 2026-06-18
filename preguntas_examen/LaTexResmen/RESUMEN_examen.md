# Resumen para el examen — Estadística geométrica en $L^2$

> **Documento de estudio consolidado** a partir de las 8 preguntas del examen, resuelto
> **inciso por inciso** y con el nivel de detalle de las soluciones largas.
>
> Todas las preguntas comparten **un mismo marco teórico** (enfoque "San Martín"):
> una base de datos es un **espacio empírico de probabilidad**; las variables son
> **funciones** sobre las etiquetas; la **esperanza condicional es la proyección
> ortogonal** en el espacio de Hilbert $L^2$. Si dominas el §0, las 8 preguntas se
> resuelven casi mecánicamente.
>
> **Cómo usar este resumen:** lee primero el §0 (la "llave"). Cada pregunta está
> dividida en sus incisos tal como aparecen en el enunciado; los recuadros `>` marcan
> las ideas que hay que poder reproducir de memoria. Al final hay tres tablas de
> repaso rápido (fórmulas, símbolos, teoremas).

---

## Índice

- [§0. Marco teórico común (la llave de todo)](#0-marco-teórico-común-la-llave-de-todo)
- [§1. Base de datos (reformas fiscales)](#1-base-de-datos-reformas-fiscales)
- [§2. Coherencia en una base de datos (falacia ecológica)](#2-coherencia-en-una-base-de-datos-falacia-ecológica)
- [§3. Modelos lineales encajonados (nested)](#3-modelos-lineales-encajonados-nested)
- [§4. Teorema de Pitágoras ($R^2$ vs $r_n^2$)](#4-teorema-de-pitágoras-r2-vs-r_n2)
- [§5. Esperanza condicional (construcción de la proyección + ley total)](#5-esperanza-condicional-construcción-de-la-proyección--ley-total)
- [§6. Igualdad del paralelogramo](#6-igualdad-del-paralelogramo)
- [§7. El término de error — dos miradas](#7-el-término-de-error--dos-miradas)
- [§8. Independencia de funciones ($n=13$ primo)](#8-independencia-de-funciones-n13-primo)
- [Apéndice — Chuleta de conceptos transversales](#apéndice--chuleta-de-conceptos-transversales)
- [Tablas de repaso rápido](#tablas-de-repaso-rápido)

---

## 0. Marco teórico común (la llave de todo)

### 0.1 Espacio empírico de probabilidad

Una base de datos es la terna $(\Omega_n, \mathcal{P}(\Omega_n), P_n)$ donde:

- $\Omega_n = \{\omega_1,\dots,\omega_n\}$ es el **conjunto de etiquetas (labels)** = los
  **individuos estadísticos**, lo que identifica a cada fila de forma única. **No** confundir
  la etiqueta con el valor de una variable: la etiqueta *identifica*, la variable *describe*.
- $\mathcal{P}(\Omega_n)$ es el $\sigma$-álgebra trivial (todas las partes): como $\Omega_n$ es
  finito, todo subconjunto es un evento medible.
- Probabilidad **uniforme** (regla de Laplace): cada etiqueta pesa $p_\omega = 1/n$, y
  $P_n(A) = |A|/n$ para todo $A\subseteq\Omega_n$ ("casos favorables sobre casos totales").
- Una **variable** es una función $X:\Omega_n \to \mathrm{Im}(X)$. **Cada columna** de la tabla
  es una variable.

### 0.2 Partición inducida por una variable

Toda variable $X$ con $\mathrm{Im}(X)=\{x_1,\dots,x_r\}$ induce una **partición** de $\Omega_n$
en sus **fibras / conjuntos de nivel / bloques**:
$$A_i := X^{-1}(\{x_i\}) = \{\omega : X(\omega)=x_i\}.$$

Que esto sea una partición se verifica en tres pasos (útil tenerlo a mano para la P8):

1. **No vacíos:** como $x_i\in\mathrm{Im}(X)$, existe $\omega$ con $X(\omega)=x_i$, luego $A_i\neq\varnothing$.
2. **Disjuntos:** $X$ es función, a cada $\omega$ le toca un único valor $\Rightarrow A_i\cap A_k=\varnothing$ si $i\neq k$.
3. **Cubren $\Omega_n$:** todo $\omega$ tiene un valor $X(\omega)$, que es alguno de los $x_i$.

La notación $\sigma(X)$ designa el **subespacio de $L^2$ de funciones que son constantes sobre
cada bloque**. Su dimensión es el número de bloques.

### 0.3 $L^2$ como espacio de Hilbert

$L^2(\Omega_n,\mathcal{P}(\Omega_n),P_n)$ con producto interno y norma:
$$\langle U,V\rangle = E[UV] = \sum_{\omega} U(\omega)V(\omega)\,P_n(\{\omega\})
= \tfrac1n\sum_{i=1}^n U(\omega_i)V(\omega_i),\qquad \|U\|^2=\langle U,U\rangle=E[U^2].$$

Como $\Omega_n$ es finito, **toda** función real está en $L^2$ y el espacio es completo
$\Rightarrow$ es un **espacio de Hilbert**. La completitud es lo que permite invocar el Teorema
de la Proyección (ver §6).

### 0.4 Esperanza condicional = proyección ortogonal

$$\boxed{\;E(Y\mid\sigma(X)) = P_{\sigma(X)}\,Y = \text{proyección ortogonal de } Y \text{ sobre } \sigma(X)\;}$$

Es el único $\hat Y\in\sigma(X)$ que minimiza $\|Y-\hat Y\|^2$ (**mejor predicción en media
cuadrática**). En cada individuo vale **el promedio de $Y$ sobre su bloque**:
$$E(Y\mid\sigma(X))(\omega)=\frac{1}{|B|}\sum_{\omega'\in B} Y(\omega'),\quad B=\text{bloque de }\omega.$$

> **Teorema de la Proyección (hecho central, se usa en TODAS las preguntas):** el residuo
> $Z=Y-\hat Y$ es ortogonal a **todo** el subespacio: $\langle Z,W\rangle=0$ para todo $W\in\sigma(X)$.

> **Dato clave que activa casi todas las identidades:** la función constante $\mathbf 1$ siempre
> está en $\sigma(X)$ (es constante en cada bloque). Como $E[Z]=\langle Z,\mathbf 1\rangle$,
> la ortogonalidad del residuo se traduce inmediatamente en $E[Z]=0$.

### 0.5 Tres consecuencias inmediatas (las usarás una y otra vez)

1. **Insesgamiento / esperanza total:** como $\mathbf 1\in\sigma(X)$ y $Z\perp\sigma(X)$:
   $$E[Z]=\langle Z,\mathbf 1\rangle=0 \;\Rightarrow\; E[E(Y\mid\sigma(X))]=E[Y].$$
2. **Pitágoras (norma):** $\|Y\|^2=\|\hat Y\|^2+\|Z\|^2$.
3. **Descomposición de varianza (ANOVA):**
   $$\mathrm{Var}(Y)=\underbrace{\mathrm{Var}(E(Y\mid\sigma(X)))}_{\text{entre bloques (SC}_\text{inter})}
   +\underbrace{E(\mathrm{Var}(Y\mid\sigma(X)))}_{\text{dentro de bloques (SC}_\text{intra})}.$$

### 0.6 Fórmulas de proyección útiles

- Sobre una recta generada por $v\neq0$: $\;P_{\sigma(v)}Y=\dfrac{\langle Y,v\rangle}{\langle v,v\rangle}\,v$.
- **Suma de proyecciones SOLO si hay ortogonalidad:** si $\mathcal M\perp\mathcal N$ entonces
  $P_{\mathcal M\oplus\mathcal N}=P_{\mathcal M}+P_{\mathcal N}$. Si **no** son ortogonales, esto
  **falla** y hay que resolver las **ecuaciones normales** (ver §3 y §7).

> **Mapa mental del marco.** Variable ⟶ partición ⟶ subespacio $\sigma(X)$ ⟶ proyección $=$
> esperanza condicional $=$ promedio por bloque. De ahí salen: residuo ortogonal, esperanza
> total, Pitágoras, ANOVA, $R^2$. Toda la asignatura es una sola idea geométrica repetida.

---

## 1. Base de datos (reformas fiscales)

**Planteo.** Un centro ligado a la Comisión Europea entrega una tabla País × Año para medir el
impacto de reformas fiscales. La tabla tiene **20 individuos** (4 países × 5 años) y 9 variables.
Marco: $\mathbb R^{20}$, pesos uniformes $1/20$, $\langle X,Y\rangle=\tfrac{1}{20}\sum_i X_iY_i$.

> **Trampa del enunciado:** menciona 5 países (incluida Italia), pero la Tabla 1 solo trae 4
> (FRA, DEU, BEL, ESP). El cardinal real de etiquetas es $20$, no $25$.

### Inciso 1 — Conjunto de etiquetas (labels)

Las etiquetas son los **identificadores de individuos**. Cada fila queda identificada por el par
**(País, Año)**, de modo que
$$I=\{(\mathrm{FRA},2015),(\mathrm{FRA},2016),\dots,(\mathrm{ESP},2019)\},\qquad |I|=4\times5=20.$$
Cada label $i\in I$ pesa $p_i=1/20$, y cada columna es una función $X:I\to\mathbb R$.
Ni el país solo ni el año solo identifican un individuo: hace falta **el par**.

### Inciso 2 — Partición inducida por FR (reforma fiscal)

FR es **binaria**, así que induce **2 clases**:
$$I=A_0\sqcup A_1,\quad A_0=\{FR=0\}\ (\text{sin reforma}),\quad A_1=\{FR=1\}\ (\text{con reforma}).$$
FR$=0$ corresponde a 2015–2016 y FR$=1$ a 2017–2019, en los 4 países:
$$|A_0|=4\times2=8,\qquad |A_1|=4\times3=12.$$

### Inciso 3 — Esperanza condicional $E(\text{Inv}\mid\text{Año})$

La variable **Año** induce **5 clases** (una por año), cada una con 4 individuos. La esperanza
condicional vale, en cada individuo, **el promedio de su clase**:
$$E(\text{Inv}\mid\text{Año})(i)=\frac{1}{|C_{\text{Año}(i)}|}\sum_{j:\ \text{Año}(j)=\text{Año}(i)}\text{Inv}(j).$$

| Año | Valores (FRA, DEU, BEL, ESP) | Suma | Promedio |
|-----|------------------------------|------|----------|
| 2015 | 19.3, 18.5, 18.0, 16.9 | 72.7 | **18.175** |
| 2016 | 17.1, 17.0, 17.1, 18.0 | 69.2 | **17.300** |
| 2017 | 15.6, 19.7, 20.2, 19.8 | 75.3 | **18.825** |
| 2018 | 21.9, 20.5, 21.2, 20.4 | 84.0 | **21.000** |
| 2019 | 16.7, 18.4, 19.4, 18.8 | 73.3 | **18.325** |

$E(\text{Inv}\mid\text{Año})$ es **constante por año** con esos valores; la media global se conserva
(insesgamiento, §0.5): $\bar x = 374.5/20 = \mathbf{18.725}$.

### Inciso 4 — Cálculo de $r^2$ e interpretación geométrica

$r^2$ es la **razón de correlación** $\eta^2$ = proporción de varianza explicada por el Año.
Por la descomposición ANOVA / Pitágoras $\mathrm{Var}(\text{Inv})=\mathrm{Var}(E(\text{Inv}\mid\text{Año}))+E(\mathrm{Var}(\text{Inv}\mid\text{Año}))$:

$$SC_\text{inter}=\sum_a n_a(\bar x_a-\bar x)^2 = 4\times7.67875 = 30.715,$$
$$SC_\text{total}=\sum_i x_i^2-20\bar x^2 = 7066.41-7012.5125 = 53.8975,$$
$$SC_\text{intra}=SC_\text{total}-SC_\text{inter}=23.1825.$$
$$\boxed{\,r^2=\frac{SC_\text{inter}}{SC_\text{total}}=\frac{30.715}{53.8975}\approx 0.570\ (\approx 57\%).\,}$$

**Interpretación geométrica.** Sea $F\subseteq\mathbb R^{20}$ el subespacio de variables
constantes por año. Entonces $E(\text{Inv}\mid\text{Año})=P_F(\text{Inv})$ es la **proyección
ortogonal** de Inv sobre $F$, y el residuo es ortogonal a $F$. Con vectores **centrados** en $\bar x$:
$$\|\text{Inv}-\bar x\|^2=\underbrace{\|\hat Y-\bar x\|^2}_{\text{explicada}}+\underbrace{\|\text{Inv}-\hat Y\|^2}_{\text{residual}}\quad(\text{Pitágoras}),$$
y $r^2=\cos^2\theta$, donde $\theta$ es el ángulo entre la variable centrada Inv y el subespacio $F$.
Aquí $\cos\theta=\sqrt{0.57}\approx0.755\Rightarrow\theta\approx41^\circ$. **El año explica ~57 %**
de la variabilidad de la inversión; el resto vive dentro de los años (entre países).

---

## 2. Coherencia en una base de datos (falacia ecológica)

**Planteo.** Datos electorales: dos candidatos; cada voto, con **folio único**, toma valor
candidato 1 / candidato 2 / nulo / blanco. Cada folio está en una mesa, que pertenece a un local,
en un distrito. Se argumenta (mal) que se pueden usar datos agregados para explicar individuos.

### Inciso 1 — ¿Cuáles son las etiquetas (labels)?

Las etiquetas son **los folios** (un voto = un individuo):
$$\Omega_N=\{\omega_1,\dots,\omega_N\},\quad N=\#\{\text{votos emitidos}\},\quad p_\omega=1/N.$$
La **opción de voto NO es una etiqueta**: es una **variable**
$$V:\Omega_N\to\{\text{cand.1, cand.2, nulo, blanco}\}$$
que a cada folio le asigna su contenido. **Confundir la etiqueta (folio) con el valor de una
variable (la opción) es el error conceptual de fondo:** el candidato no individualiza la unidad,
solo la describe. Para ser etiquetas legítimas, los folios deben ser **únicos, exhaustivos y
mutuamente excluyentes**.

### Inciso 2 — Mesa, Local y Distrito, ¿son variables? ¿Qué particiones inducen?

**Sí:** las tres son funciones bien definidas sobre los folios:
$$M:\Omega_N\to\{\text{mesas}\},\quad L:\Omega_N\to\{\text{locales}\},\quad D:\Omega_N\to\{\text{distritos}\}.$$
Cada una induce una partición por preimágenes: $M^{-1}(m)$ = folios de la mesa $m$, etc.

Lo esencial es que estas particiones son **jerárquicas / encajonadas**: cada mesa está en un
único local y cada local en un único distrito. Por tanto la partición de Mesa **refina** a la de
Local, que refina a la de Distrito. En términos de subespacios:
$$\boxed{\;\sigma(D)\subseteq\sigma(L)\subseteq\sigma(M)\subseteq L^2(\Omega_N).\;}$$

| Variable | Bloques (preimágenes) | Finura relativa |
|----------|-----------------------|-----------------|
| Mesa | folios de una misma mesa | **más fina** |
| Local | unión de mesas del local | intermedia |
| Distrito | unión de locales del distrito | **más gruesa** |

Equivalentemente: todo bloque de Distrito es unión de bloques de Local, y todo bloque de Local
es unión de bloques de Mesa. Existe la función de paso $\text{mesa}\mapsto\text{local}\mapsto\text{distrito}$,
**pero no al revés** (más fino determina más grueso, no a la inversa).

### Inciso 3 — Condiciones para incorporar la tasa de cesantía comunal $C$

La tasa de cesantía es un **dato agregado a nivel de comuna**, no de folio. Para ser una variable
de esta base, $C$ debe poder escribirse como función sobre las etiquetas,
$C:\Omega_N\to\mathbb R$, $\omega\mapsto$ tasa de cesantía de la comuna del folio $\omega$.
Esto exige **dos condiciones**:

1. **Coherencia con la jerarquía (que la comuna sea variable bien definida):** debe existir una
   asignación unívoca folio → comuna, lo cual ocurre solo si cada local (o mesa/distrito)
   pertenece a una única comuna. Así queda definida $K:\Omega_N\to\{\text{comunas}\}$ compatible
   con $\sigma(D)\subseteq\sigma(L)\subseteq\sigma(M)$. Si un local cayera en dos comunas, $K$ no
   estaría bien definida.
2. **Medibilidad:** $C$ debe ser **constante por comuna**, es decir
   $$\boxed{\,C\in\sigma(K),\quad K\text{ coherente con la jerarquía territorial.}\,}$$

### El punto de fondo — ¿por qué "falacia ecológica"?

Como $C\in\sigma(K)$ es **constante dentro de cada comuna**, solo aporta variabilidad **entre**
comunas, no dentro. Toda regresión de $V$ sobre $C$ **factoriza** por la proyección comunal:
$$E(V\mid C)=E(V\mid\sigma(K)),$$
y no distingue a dos votantes de la misma comuna. La descomposición lo deja explícito:
$$\mathrm{Var}(V)=\underbrace{\mathrm{Var}(E(V\mid\sigma(K)))}_{\text{entre comunas (lo que ``ve'' }C)}
+\underbrace{E(\mathrm{Var}(V\mid\sigma(K)))}_{\text{dentro: invisible para }C}.$$

> **Falacia ecológica (Robinson, 1950):** inferir el comportamiento individual a partir de
> correlaciones medidas a nivel agregado. Una asociación fuerte *entre* comunas no dice nada de
> lo que ocurre *dentro*, que es justo el nivel individual de interés. La relación incluso puede
> **invertirse** (**paradoja de Simpson**).

**Ejemplo de inversión (Simpson).** Dos comunas; a nivel comunal "más cesantía → más voto al
candidato 1", pero dentro de cada comuna los empleados votan **menos** al candidato 1:

| Comuna | Cesantía (agregada) | % voto cand. 1 (individual) |
|--------|--------------------|------------------------------|
| A — cesantes | 20% | 70% |
| A — empleados | 20% | 40% |
| B — cesantes | 5% | 55% |
| B — empleados | 5% | 25% |

**Moraleja.** Incorporar $C$ es legítimo y útil para describir diferencias *entre* comunas, pero
**no es información a nivel de folio**. Para una explicación individual hace falta una variable
definida realmente sobre el folio (p. ej. la condición de empleo de cada votante) que induzca una
partición **más fina que la comunal**.

---

## 3. Modelos lineales encajonados (nested)

**Planteo.** Modelos $M_1: E(Y\mid X_1)=\alpha_0+\alpha_1X_1$ y
$M_2: E(Y\mid X_1,X_2)=\alpha_0+\alpha_1X_1+\alpha_2X_2$. Se dice que $M_1$ está encajonado en
$M_2$ (poniendo $\alpha_2=0$), lo que justifica la práctica de "ajustar $M_1$, agregar $X_2$ y ver
cuál ajusta mejor". El ejercicio pide **dos ejemplos** y una **digresión crítica**.

**Marco.** $\Omega_3=\{\omega_1,\omega_2,\omega_3\}$, uniforme $1/3$, vectores en $\mathbb R^3$,
$\langle U,V\rangle=\tfrac13\sum_i U_iV_i$. Hecho clave a justificar (no darlo por hecho):
$P_{\mathcal M\oplus\mathcal N}=P_{\mathcal M}+P_{\mathcal N}$ **solo si** $\mathcal M\perp\mathcal N$.
Y la fórmula de proyección sobre una recta: $P_{\sigma(v)}Y=\frac{\langle Y,v\rangle}{\langle v,v\rangle}v$.

### Ejemplo 1 — regresores ORTOGONALES (se cumple la aditividad)

Vectores: $\mathbf 1=(1,1,1)$, $X=(-1,0,1)$, $Y=(y_1,y_2,y_3)$.

**Paso 1 — productos internos.** $\langle\mathbf 1,\mathbf 1\rangle=1$,
$\langle X,X\rangle=\tfrac23$, $\langle\mathbf 1,X\rangle=\tfrac13(-1+0+1)=0$,
$\langle Y,\mathbf 1\rangle=\bar y$, $\langle Y,X\rangle=\tfrac{y_3-y_1}{3}$.

**Paso 2 — observación central:** $\langle\mathbf 1,X\rangle=0\Rightarrow \mathbf 1\perp X$.

**Paso 3 — proyecciones simples.** $P_{\sigma(\mathbf 1)}Y=\bar y\,\mathbf 1$ y
$P_{\sigma(X)}Y=\frac{(y_3-y_1)/3}{2/3}X=\frac{y_3-y_1}{2}X$.

**Paso 4 — proyección sobre el plano.** Como $\{\mathbf 1,X\}$ es base ortogonal de $\sigma(\mathbf 1,X)$:
$$E(Y\mid\sigma(\mathbf 1,X))=\bar y\,\mathbf 1+\tfrac{y_3-y_1}{2}X
= E(Y\mid\sigma(\mathbf 1))+E(Y\mid\sigma(X)).\;\checkmark$$
Por componentes: $E(Y\mid\sigma(\mathbf 1,X))=\big(\bar y-\tfrac{y_3-y_1}{2},\ \bar y,\ \bar y+\tfrac{y_3-y_1}{2}\big)$
— la recta de regresión por MCO evaluada en cada punto.

> **¿A qué se debe la igualdad?** A la **ortogonalidad** $\langle\mathbf 1,X\rangle=0$. Cuando los
> ejes son ortogonales, proyectar sobre el plano = sumar las proyecciones sobre cada eje, sin
> interferencia. No hay magia: es Pitágoras / descomposición ortogonal en acción.

### Ejemplo 2 — regresores NO ortogonales (falla la aditividad)

Ahora $\mathbf 1=(1,1,1)$, $Z=(1,0,1)$, $Y=(y_1,y_2,y_3)$.

**Paso 1 — productos internos.** $\langle\mathbf 1,\mathbf 1\rangle=1$, $\langle Z,Z\rangle=\tfrac23$,
$\langle\mathbf 1,Z\rangle=\tfrac23$, $\langle Y,\mathbf 1\rangle=\bar y$, $\langle Y,Z\rangle=\tfrac{y_1+y_3}{3}$.

**Paso 2 — observación central:** $\langle\mathbf 1,Z\rangle=\tfrac23\neq0\Rightarrow$ **no** ortogonales.
No estamos autorizados a sumar proyecciones; hay que resolver **ecuaciones normales**.

**Paso 3 — ecuaciones normales.** Buscamos $a\mathbf 1+bZ$ con residuo $\perp\mathbf 1$ y $\perp Z$:
$$\begin{cases} a\langle\mathbf 1,\mathbf 1\rangle+b\langle Z,\mathbf 1\rangle=\langle Y,\mathbf 1\rangle\\ a\langle\mathbf 1,Z\rangle+b\langle Z,Z\rangle=\langle Y,Z\rangle\end{cases}
\Rightarrow\begin{cases} a+\tfrac23 b=\tfrac{y_1+y_2+y_3}{3}\\ \tfrac23 a+\tfrac23 b=\tfrac{y_1+y_3}{3}\end{cases}$$
Restando: $(1-\tfrac23)a=\tfrac{y_2}{3}\Rightarrow a=y_2$; luego $b=\tfrac{y_1+y_3}{2}-y_2$. Por tanto
$$E(Y\mid\sigma(\mathbf 1,Z))=\Big(\tfrac{y_1+y_3}{2},\,y_2,\,\tfrac{y_1+y_3}{2}\Big).$$

**Paso 4 — suma de proyecciones simples.**
$E(Y\mid\sigma(\mathbf 1))=\bar y\,\mathbf 1$ y $E(Y\mid\sigma(Z))=\frac{(y_1+y_3)/3}{2/3}Z=\frac{y_1+y_3}{2}Z$.

**Paso 5 — comparación (mirar $\omega_2$):**
$$[E(Y\mid\sigma(\mathbf 1,Z))]_{\omega_2}=y_2 \quad\text{vs.}\quad [E(Y\mid\sigma(\mathbf 1))+E(Y\mid\sigma(Z))]_{\omega_2}=\tfrac{y_1+y_2+y_3}{3}.$$
Coinciden **solo si** $2y_2=y_1+y_3$ (caso particular). En general:
$$E(Y\mid\sigma(\mathbf 1,Z))\neq E(Y\mid\sigma(\mathbf 1))+E(Y\mid\sigma(Z)).$$
Al no ser ortogonales, $\sigma(\mathbf 1)$ y $\sigma(Z)$ se **solapan**; sumar proyecciones cuenta
dos veces la parte común. El coeficiente correcto $b$ corrige ese solapamiento.

### Pregunta — ¿cuándo son genuinamente "encajonados" $M_1$ y $M_2$?

Hay que distinguir **dos nociones** que la práctica confunde:

1. **Encajonamiento como inclusión de subespacios.** Siempre se cumple
   $\sigma(\mathbf 1)\subseteq\sigma(\mathbf 1,X_1)\subseteq\sigma(\mathbf 1,X_1,X_2)$, sin condición.
   Esto valida **solo una cosa:** comparar ajustes ($M_2$ nunca ajusta peor, pues proyecta sobre
   un espacio más grande: $\|Y-P_{\sigma(\mathbf 1,X_1)}Y\|^2-\|Y-P_{\sigma(\mathbf 1,X_1,X_2)}Y\|^2\ge0$).
2. **Encajonamiento como descomposición aditiva de la contribución de cada predictor** (lo que la
   práctica realmente usa: "agrego $X_2$ y le atribuyo la mejora sin que cambie lo demás"). Esto
   se cumple **si y solo si** los regresores son ortogonales:
   $$E(Y\mid\sigma(\mathbf 1,X))=E(Y\mid\sigma(\mathbf 1))+E(Y\mid\sigma(X))\iff \langle\mathbf 1,X\rangle=0.$$

> **Criterio crítico.** $M_1$ y $M_2$ son "encajonados" en **sentido fuerte** (el coeficiente de
> $X_1$ no cambia al añadir $X_2$ y la mejora se atribuye limpiamente a $X_2$) **solo cuando
> $X_2\perp\sigma(\mathbf 1,X_1)$**. Si los predictores están correlacionados, agregar $X_2$
> **altera** el coeficiente de $X_1$. La forma correcta de "encajonar" es **ortogonalizar primero**:
> comparar $Y$ contra la parte de $X_2$ residual respecto de $\sigma(\mathbf 1,X_1)$ — esta es la
> idea del **teorema de Frisch–Waugh–Lovell**. La receta ingenua "ajusto, agrego, comparo" solo es
> interpretable bajo ortogonalidad; en general puede ser engañosa.

---

## 4. Teorema de Pitágoras ($R^2$ vs $r_n^2$)

**Setup.** $X,Y\in L^2(\Omega_n,\mathcal P(\Omega_n),P_n)$; $\sigma(X)$ = funciones constantes por
bloque. Con $\hat Y=E(Y\mid\sigma(X))$, $Z=Y-\hat Y$, el Teorema de la Proyección da la
descomposición $Y=\hat Y+Z$ con $Z\perp\sigma(X)$, de donde (1) $\|Y\|^2=\|\hat Y\|^2+\|Z\|^2$.
Recordatorio: $\mathbf 1\in\sigma(X)\Rightarrow\langle Z,\mathbf 1\rangle=E[Z]$.

### Inciso 1 — Demostrar $E(Y-E(Y\mid\sigma(X)))=0$

Queremos $E[Z]=0$. **Paso 1:** $\mathbf 1\in\sigma(X)$. **Paso 2:** Teorema de la Proyección con
$W=\mathbf 1$: $\langle Z,\mathbf 1\rangle=0$. **Paso 3:** por definición del producto interno,
$\langle Z,\mathbf 1\rangle=E[Z\cdot\mathbf 1]=E[Z]$. **Paso 4:** combinando, $E[Z]=0$. $\blacksquare$

**Consecuencia (para el inciso 2):** como $Z=Y-\hat Y$ y la esperanza es lineal,
$0=E[Z]=E[Y]-E[\hat Y]\Rightarrow E[\hat Y]=E[Y]=:\mu$ (insesgamiento de la esperanza condicional).

### Inciso 2 — Probar que (1) equivale a la descomposición de varianzas (2)

(2) es $\mathrm{Var}(Y)=\mathrm{Var}[E(Y\mid\sigma(X))]+\mathrm{Var}[Y-E(Y\mid\sigma(X))]$.
Identidad básica: $\mathrm{Var}(U)=E[U^2]-[E(U)]^2$.

**Paso 1 — esperanzas:** del inciso 1, $E[Y]=\mu$, $E[\hat Y]=\mu$, $E[Z]=0$.

**Paso 2 — varianzas:** $\mathrm{Var}(Y)=\|Y\|^2-\mu^2$, $\mathrm{Var}(\hat Y)=\|\hat Y\|^2-\mu^2$,
$\mathrm{Var}(Z)=\|Z\|^2$ (porque $E[Z]=0$).

**Paso 3 — de (1) a (2):** partimos de $E[Y^2]=E[\hat Y^2]+E[Z^2]$ y **restamos $\mu^2$** a ambos
lados (lícito, misma cantidad). Agrupamos el $-\mu^2$ con el **primer** término (es $\hat Y$ quien
lleva la media, $Z$ no):
$$\underbrace{E[Y^2]-\mu^2}_{\mathrm{Var}(Y)}=\underbrace{(E[\hat Y^2]-\mu^2)}_{\mathrm{Var}(\hat Y)}+\underbrace{E[Z^2]}_{\mathrm{Var}(Z)}.$$

**Paso 4 — de (2) a (1):** "restar $\mu^2$" es reversible; sumando $\mu^2$ recuperamos (1). Como
cada una se obtiene de la otra por una operación reversible, **(1) ⟺ (2)**. $\blacksquare$
La equivalencia descansa íntegramente en el inciso 1 ($E[\hat Y]=E[Y]$, así que el ajuste $-\mu^2$
es común a ambos lados).

### Inciso 3 — Comparar $R^2$ y $r_n^2$

Con $a:=\|\hat Y\|^2=E[\hat Y^2]$ y $b:=\|Y\|^2=E[Y^2]$:
$$r_n^2=\frac{\|E(Y\mid\sigma(X))\|^2}{\|Y\|^2}=\frac{a}{b}\quad\text{(momentos crudos)},\qquad
R^2=\frac{\mathrm{Var}(\hat Y)}{\mathrm{Var}(Y)}=\frac{a-\mu^2}{b-\mu^2}\quad\text{(varianzas)}.$$
La única diferencia: $R^2$ resta $\mu^2$ arriba y abajo (segundos momentos **centrados**).

**Ordenamiento de $\mu^2, a, b$:** de (1), $b=a+E[Z^2]\ge a$ ($\Rightarrow a\le b$); por Jensen
(o $\mathrm{Var}(\hat Y)\ge0$), $a=E[\hat Y^2]\ge[E\hat Y]^2=\mu^2$. Luego $\mu^2\le a\le b$
(suponemos $\mathrm{Var}(Y)>0$, i.e. $b>\mu^2$, para que $R^2$ esté definido).

**Comparación por producto cruzado** (denominadores positivos):
$$R^2\le r_n^2\iff (a-\mu^2)b\le a(b-\mu^2)\iff ab-\mu^2 b\le ab-\mu^2 a\iff \mu^2 a\le\mu^2 b,$$
cierto porque $\mu^2\ge0$ y $a\le b$. Por tanto
$$\boxed{\,R^2\le r_n^2.\,}$$
**Igualdad** sii $\mu^2(b-a)=0$: o bien $\mu=E[Y]=0$ (variable **centrada**), o bien $a=b$
($Z=0$, i.e. $Y=\hat Y$ c.s., **ajuste perfecto**).

**Interpretación.** $R^2$ es el **genuino coeficiente de determinación** (proporción de varianza
explicada, **invariante a traslaciones**); $r_n^2$ usa momentos crudos, incorpora la media y
**sobrestima** el ajuste. La diferencia crece con $|\mu|$.

**Ejemplo concreto.** $\Omega=\{1,2,3,4\}$ uniforme, bloques $\{1,2\},\{3,4\}$, $Y=(0,2,2,4)$.
Entonces $\hat Y=(1,1,3,3)$ (media por bloque), $\mu=2$, $b=E[Y^2]=6$, $a=E[\hat Y^2]=5$,
$\mathrm{Var}(Y)=2$, $\mathrm{Var}(\hat Y)=1$, $Z=(-1,1,-1,1)$ con $E[Z]=0$, $\mathrm{Var}(Z)=1$
(verifica (2): $2=1+1$ ✓). Resultado:
$$r_n^2=\tfrac56\approx0.833,\qquad R^2=\tfrac12=0.5\quad(R^2<r_n^2\ \checkmark).$$
$r_n^2$ sugiere un "ajuste" mucho mayor solo porque no descuenta la media $\mu=2$.

---

## 5. Esperanza condicional (construcción de la proyección + ley total)

**Setup.** $(H,\langle\cdot,\cdot\rangle)$ Hilbert, $M\subset H$ subespacio **cerrado**. Dato de
partida: para cada $x\in H$ existe un único $y_0\in M$ con $\|x-y_0\|=\inf_{y\in M}\|x-y\|$.
Herramienta repetida — **identidad de Pitágoras:** si $u\perp v$, $\|u+v\|^2=\|u\|^2+\|v\|^2$
(el cruzado $2\langle u,v\rangle$ desaparece *solo* por la ortogonalidad).

### Inciso 1 — Construir $P:H\to M$ y probar que es proyección ortogonal

**Definición.** $Px:=y_0$, el minimizador (bien definido por la **unicidad** del dato).

**Lema clave (caracterización por ortogonalidad).** $\;y_0$ minimiza $\iff (x-y_0)\perp M$.
- ($\Leftarrow$) Si $x-y_0\perp M$, para todo $y\in M$ se tiene $y_0-y\in M$, y por Pitágoras
  $\|x-y\|^2=\|x-y_0\|^2+\|y_0-y\|^2\ge\|x-y_0\|^2$, con igualdad sii $y=y_0$.
- ($\Rightarrow$) Si $y_0$ minimiza, fijo $m\in M$ y considero
  $\varphi(t)=\|x-(y_0+tm)\|^2=\|x-y_0\|^2-2t\langle x-y_0,m\rangle+t^2\|m\|^2$. Tiene mínimo en
  $t=0\Rightarrow\varphi'(0)=-2\langle x-y_0,m\rangle=0$. Como $m$ es arbitrario, $x-y_0\perp M$.

A partir del lema salen casi solas las propiedades:
- **Lineal:** $(ax_1+bx_2)-(aPx_1+bPx_2)\perp M\Rightarrow$ por unicidad $P(ax_1+bx_2)=aPx_1+bPx_2$.
- **Idempotente:** $x\in M\Rightarrow Px=x$ (distancia 0), luego $P^2=P$.
- **Autoadjunto:** $\langle Px,z\rangle=\langle Px,Pz\rangle=\langle x,Pz\rangle$ (usando
  $z-Pz\perp M\ni Px$ y $x-Px\perp M\ni Pz$). Junto con $P^2=P$, esto **define** proyección ortogonal.
- **Contracción:** $\|x\|^2=\|Px\|^2+\|x-Px\|^2\ge\|Px\|^2\Rightarrow\|P\|\le1$.
- **Descomposición:** $H=M\oplus M^\perp$ (núcleo de $P$ es $M^\perp$).

### Inciso 2 — ¿Por qué la proyección ortogonal ES la esperanza condicional?

Reconocemos el Hilbert probabilístico $H=L^2(\Omega,\mathcal F,P)$ con $\langle X,Y\rangle=E[XY]$
(la norma $\|X\|^2=E[X^2]$ mide el **error cuadrático medio**). Dada una sub-$\sigma$-álgebra
$\mathcal G\subseteq\mathcal F$, el subespacio $M=L^2(\Omega,\mathcal G,P)$ ($\mathcal G$-medibles)
es lineal y **cerrado**. La esperanza condicional $E[X\mid\mathcal G]$ se define como la variable
$Z\in M$ ($\mathcal G$-medible) con
$$E[(X-Z)W]=0\quad\forall W\in M.$$
Pero $E[(X-Z)W]=\langle X-Z,W\rangle$, así que la condición dice **exactamente** $X-Z\perp M$. Por
el lema del inciso 1, esto caracteriza $Z=PX$. Por tanto
$$\boxed{\,E[X\mid\mathcal G]=PX=\text{proyección ortogonal de }X\text{ sobre }M.\,}$$
Interpretación: la **mejor predicción de $X$ en media cuadrática** usando solo la información de $\mathcal G$.

### Inciso 3 — Demostrar $E[E(Y\mid\sigma(X))]=E(Y)$ (ley de la esperanza total)

Concretamos: $\Omega_n$ finito, $H=L^2(\Omega_n)$, $X$ induce bloques $B_1,\dots,B_k$, y
$M=\sigma(X)$ (constantes por bloque), cerrado por dimensión finita.

**Observación clave:** la constante $\mathbf 1\in M$, y para toda $Z\in H$,
$E[Z]=\sum_\omega Z(\omega)P_n(\{\omega\})=\langle Z,\mathbf 1\rangle$.

**Aplicar Proyección:** $\hat Y:=E[Y\mid\sigma(X)]$ es la proyección de $Y$ sobre $M$, luego
$Y-\hat Y\perp M$. **Conclusión:** como $\mathbf 1\in M$,
$$E[Y]-E[\hat Y]=\langle Y,\mathbf 1\rangle-\langle\hat Y,\mathbf 1\rangle=\langle Y-\hat Y,\mathbf 1\rangle=0\;\Rightarrow\;E[\hat Y]=E[Y].\;\blacksquare$$

**Verificación explícita (opcional).** Como $Y-\hat Y\perp M$, tomando $\mathbf 1_{B_j}\in M$ se
obtiene $\hat Y|_{B_j}=\frac{1}{P_n(B_j)}\sum_{\omega\in B_j}Y(\omega)P_n(\{\omega\})=E[Y\mid B_j]$.
Entonces $E[\hat Y]=\sum_j P_n(B_j)E[Y\mid B_j]=\sum_j\sum_{\omega\in B_j}Y(\omega)P_n(\{\omega\})=E[Y]$,
que es la ley de la esperanza total.

---

## 6. Igualdad del paralelogramo

**Enunciado.** En un espacio con producto interno, $\|x+y\|^2+\|x-y\|^2=2\|x\|^2+2\|y\|^2$.

### Inciso 1 — Demostración

Marco: $V$ real con producto interno y norma inducida $\|v\|^2=\langle v,v\rangle$ (**único** punto
donde se usa que la norma viene de un producto interno).

**Paso 1 — expandir $\|x+y\|^2$** por bilinealidad y simetría ($\langle x,y\rangle=\langle y,x\rangle$):
$$\|x+y\|^2=\langle x+y,x+y\rangle=\|x\|^2+2\langle x,y\rangle+\|y\|^2.$$
**Paso 2 — expandir $\|x-y\|^2$** análogamente:
$$\|x-y\|^2=\|x\|^2-2\langle x,y\rangle+\|y\|^2.$$
**Paso 3 — sumar:** los cruzados $\pm2\langle x,y\rangle$ **se cancelan**:
$$\|x+y\|^2+\|x-y\|^2=2\|x\|^2+2\|y\|^2.\;\blacksquare$$
(Caso complejo: simetría hermitiana, cruzados $=2\,\mathrm{Re}\langle x,y\rangle$, cancela igual.)

### Inciso 2 — Relevancia para definir la esperanza condicional

**Idea central:** la esperanza condicional es una **proyección ortogonal en $L^2$**, y la igualdad
del paralelogramo es la herramienta que garantiza que esa proyección **existe y es única**
(Teorema de la Proyección). Esquema de la prueba (lo que el examen quiere ver):

- $L^2(\Omega,\mathcal F,P)$ con $\langle X,Y\rangle=E[XY]$ es Hilbert (norma viene de producto
  interno ⇒ cumple la igualdad; y es completo). $M=L^2(\Omega,\mathcal G,P)$ es subespacio cerrado.
- $E[X\mid\mathcal G]$ es el $Y^*\in M$ que minimiza $\|X-Y\|^2$ (mejor predicción cuadrática).
- **Existencia del minimizador (aquí entra la igualdad):** sea $d=\inf_{Y\in M}\|X-Y\|$ y $(Y_n)\subset M$
  minimizante ($\|X-Y_n\|\to d$). Aplicando la igualdad a $a=X-Y_n$, $b=X-Y_m$:
  $$\|Y_m-Y_n\|^2=2\|X-Y_n\|^2+2\|X-Y_m\|^2-4\Big\|X-\tfrac{Y_n+Y_m}{2}\Big\|^2.$$
  Como $\tfrac{Y_n+Y_m}{2}\in M$, $\|X-\tfrac{Y_n+Y_m}{2}\|^2\ge d^2$, luego
  $$\|Y_m-Y_n\|^2\le 2\|X-Y_n\|^2+2\|X-Y_m\|^2-4d^2\xrightarrow[n,m\to\infty]{}2d^2+2d^2-4d^2=0.$$
  Entonces $(Y_n)$ es de **Cauchy**.
- **Cierre:** por completitud de $L^2$, $Y_n\to Y^*$; por ser $M$ cerrado, $Y^*\in M$; por
  continuidad de la norma, $\|X-Y^*\|=d$. Existe (y se prueba único) el minimizador $=E[X\mid\mathcal G]$.

> **Conclusión.** Sin la estructura de producto interno —es decir, sin que la norma satisfaga la
> igualdad del paralelogramo— el argumento de Cauchy se rompe y no se garantiza la mejor
> aproximación. Por eso la esperanza condicional se define en $L^2$ y **no** con la norma del
> supremo (que sirve, p. ej., para Kolmogorov–Smirnov, pero no proviene de un producto interno).

### Inciso 3 — La norma del supremo en $\mathbb R^n$ ($n\ge2$) NO la satisface

Contraejemplo: $x=e_1=(1,0,\dots)$, $y=e_2=(0,1,0,\dots)$.
$\|x\|_\infty=\|y\|_\infty=1$; $x+y=(1,1,0,\dots)$ y $x-y=(1,-1,0,\dots)$ dan $\|x+y\|_\infty=\|x-y\|_\infty=1$.
- Lado izquierdo: $1^2+1^2=2$.
- Lado derecho: $2\cdot1^2+2\cdot1^2=4$.

$2\neq4\Rightarrow$ falla $\Rightarrow$ no proviene de un producto interno. (En $n=1$ sí coincide con $|\cdot|$.)

### Inciso 4 — La norma del supremo en $C[0,1]$ NO la satisface

Contraejemplo: $F(x)=1$ (constante), $G(x)=x$, ambas continuas en $[0,1]$.
$\|F\|_\infty=\|G\|_\infty=1$; $(F+G)(x)=1+x\Rightarrow\|F+G\|_\infty=2$; $(F-G)(x)=1-x\Rightarrow\|F-G\|_\infty=1$.
- Lado izquierdo: $2^2+1^2=5$.
- Lado derecho: $2\cdot1^2+2\cdot1^2=4$.

$5\neq4\Rightarrow$ falla. El estadístico KS usa $\|\cdot\|_\infty$, pero esa norma no induce geometría de Hilbert.

---

## 7. El término de error — dos miradas

**Pregunta.** Contrastar las dos miradas del término de error —la de **Gauss** y la del **Teorema
de la Proyección**—, enfatizando el **estilo de análisis** y el **desarrollo metodológico** de cada una.

**Punto de partida común.** Modelo lineal $y=X\beta+\varepsilon$, con $y\in\mathbb R^n$,
$X\in\mathbb R^{n\times p}$, $\mathrm{rango}(X)=p$ (⇒ $X^\top X$ invertible), $\beta\in\mathbb R^p$,
$\varepsilon\in\mathbb R^n$. Ambas miradas llegan al **mismo** estimador
$\hat\beta=(X^\top X)^{-1}X^\top y$ por caminos opuestos. La clave: el **mismo** objeto formal
$\varepsilon$ recibe **dos significados distintos**.

### Mirada de Gauss (probabilística): el error como perturbación aleatoria con ley

- **Qué se postula.** $\varepsilon$ es un **objeto primitivo**: variable aleatoria que existe antes
  de estimar nada (causas no modeladas, error de medición, factores omitidos). No se deriva, **se supone**.
- **Supuestos Gauss–Markov.** $E[\varepsilon]=0$ (no sesga: toda la señal está en $X\beta$) y
  $\mathrm{Var}(\varepsilon)=\sigma^2 I_n$, que encierra **dos** hipótesis: diagonal constante =
  **homocedasticidad**, y ceros fuera de la diagonal = **no correlación**. Es una hipótesis sobre
  la **población**. Para inferir se agrega $\varepsilon\sim\mathcal N(0,\sigma^2 I_n)$.
- **Argumento histórico (por qué normal y por qué MCO).** (1) la media aritmética es el valor más
  plausible de una cantidad medida; (2) ¿qué ley hace que la media sea máximo-verosímil?; (3)
  resolviendo esa condición funcional, la densidad debe ser la **normal**; (4) maximizar la
  verosimilitud normal $\equiv$ **minimizar $\sum\varepsilon_i^2$** ⇒ justifica MCO desde la
  probabilidad. Más tarde (Theoria Combinationis), **Gauss–Markov**: sin normalidad, bajo solo
  $E[\varepsilon]=0$ y $\mathrm{Var}(\varepsilon)=\sigma^2 I$, $\hat\beta$ es **BLUE** (mínima
  varianza entre lineales insesgados).
- **Cómo se obtiene el estimador.** Minimizar $S(\beta)=\|y-X\beta\|^2$, derivar:
  $\frac{\partial S}{\partial\beta}=-2X^\top(y-X\beta)=0\Rightarrow X^\top X\hat\beta=X^\top y\Rightarrow\hat\beta=(X^\top X)^{-1}X^\top y$.
  El motor es la **ley de probabilidad**, no la geometría.
- **Estilo.** Paramétrico, inferencial, muestral. Productos típicos: insesgadez $E[\hat\beta]=\beta$,
  $\mathrm{Var}(\hat\beta)=\sigma^2(X^\top X)^{-1}$, tests $t$ y $F$, intervalos de confianza,
  Cramér–Rao, eficiencia. Mentalidad: **de la probabilidad hacia las propiedades**.

### Mirada de la Proyección (geométrica): el error como complemento ortogonal

- **Cambio de escenario.** No se parte de una ley sino de un espacio con producto interno.
  $\mathbb R^n$ con $\langle u,v\rangle=u^\top v$; subespacio de regresores $\mathcal C(X)=\{X\beta\}$.
- **Teorema de la Proyección.** Para todo $y$ existe descomposición única $y=\hat y+e$ con
  $\hat y\in\mathcal C(X)$, $e\perp\mathcal C(X)$, y $\hat y$ es el punto de $\mathcal C(X)$ más
  cercano a $y$. **El error $e$ se define como el complemento ortogonal** (la parte de $y$ que no
  es combinación de columnas de $X$): no se supone su ley, **se construye**.
- **Ecuaciones normales sin derivar.** $e\perp\mathcal C(X)\iff X^\top e=0\iff X^\top X\hat\beta=X^\top y$.
  Mismas ecuaciones que en Gauss, pero ahora como **traducción algebraica de la ortogonalidad**.
  Matriz de proyección $P=X(X^\top X)^{-1}X^\top$ ($P^2=P$, $P^\top=P$), $\hat y=Py$, $e=(I-P)y$.
  Minimizar $\|y-X\beta\|^2$ y proyectar coinciden porque la proyección **es** el punto más cercano.
- **Pitágoras = ANOVA sin distribución.** $\|y\|^2=\|\hat y\|^2+\|e\|^2\Rightarrow SCT=SCR+SCE$.
  La ANOVA **no requiere normalidad**; es pura geometría. La normalidad se agrega después, solo
  para ponerle una distribución a los cocientes (tests $F$).
- **Versión profunda en $L^2$.** En $L^2(\Omega)$ con $\langle U,V\rangle=E[UV]$, la ortogonalidad
  equivale a **no correlación** (con medias cero, $\langle U,V\rangle=\mathrm{Cov}(U,V)$).
  $E[Y\mid X]$ es la proyección de $Y$; el error $\varepsilon=Y-E[Y\mid X]$ es **automáticamente**
  ortogonal (no correlacionado) con los regresores. Lo que en Gauss era hipótesis, **aquí es
  teorema gratuito**.
- **Estilo.** Geométrico, estructural, libre de distribución. Productos típicos: matrices de
  proyección, Pitágoras/ANOVA, esperanza condicional como proyección, unificación
  regresión–ANOVA–predicción. Mentalidad: **de la estructura del espacio hacia las ecuaciones**;
  la probabilidad se incorpora después, solo para inferir.

### Contraste explícito

| Criterio | Mirada de Gauss | Mirada de la Proyección |
|----------|-----------------|--------------------------|
| Naturaleza del error | Primitivo (v.a. con ley supuesta) | Derivado (residuo ortogonal construido) |
| No correlación / ortogonalidad | **Supuesto** sobre la población | **Consecuencia** automática de proyectar |
| Origen de las ecuaciones normales | Anular $\partial S/\partial\beta$ | Traducir $X^\top e=0$ |
| Por qué MCO | Verosimilitud normal / BLUE | Proyección = punto más cercano |
| Papel de la distribución | Central, desde el inicio | Postergada; solo para inferir |
| Estilo metodológico | Inferencial, paramétrico, muestral | Geométrico, estructural, libre de distribución |
| Dirección del razonamiento | De la probabilidad a las propiedades | De la geometría a las ecuaciones |

### Síntesis

Convergen en el mismo $\hat\beta=(X^\top X)^{-1}X^\top y$, por caminos opuestos:
- **Gauss responde "por qué es óptimo"** (aporte **inferencial**): dada una ley del error, MCO es
  máximo-verosímil / mínima varianza.
- **La proyección responde "por qué tiene esa forma"** (aporte **estructural**): estimar es separar
  $y$ en lo que vive en $\mathcal C(X)$ y lo ortogonal.

> **Distinción fina (cerrar sin salto lógico):** los **errores verdaderos** $\varepsilon$ (objeto de
> Gauss) **no** son los **residuos** $e$ (objeto de la proyección); $e=(I-P)\varepsilon$ es una
> *estimación* de $\varepsilon$. La síntesis moderna usa ambas: la geometría da la forma y las
> propiedades algebraicas; la probabilidad de Gauss añade la capa de inferencia (distribuciones,
> tests, intervalos) sobre esa estructura.

---

## 8. Independencia de funciones ($n=13$ primo)

**Setup.** $\Omega_{13}$ empírico uniforme, $X:\Omega_{13}\to\{x_1,\dots,x_r\}$ y
$Y:\Omega_{13}\to\{y_1,\dots,y_q\}$. Cada variable induce una partición en fibras (no vacías,
disjuntas, que cubren), como en §0.2.

### Inciso 1 — ¿Se puede elegir $r$ de manera arbitraria?

**No:** hay una restricción que depende de $n$, a saber $\boxed{1\le r\le n}$. Como los $r$ bloques
son disjuntos y cubren $\Omega_n$, $n=\sum_{i=1}^r|A_i|$; y como cada $|A_i|\ge1$,
$n=\sum_i|A_i|\ge\sum_i 1=r$, luego $r\le n$. Además $r\ge1$ (la imagen de una función con dominio
no vacío tiene al menos un elemento). Casos extremos:
- $r=1\iff X$ **constante** (un solo bloque $A_1=\Omega_n$).
- $r=n\iff X$ **inyectiva** (cada bloque es un único elemento).
- cualquier $1<r<n$ es alcanzable. Con $n=13$: $r\in\{1,2,\dots,13\}$.

Intuición: no se pueden producir más valores distintos que etiquetas disponibles (cada valor
"consume" al menos una etiqueta).

### Inciso 2 — Partición inducida por el vector $(X,Y)$

Para el vector $(X,Y):\Omega_n\to\mathrm{Im}(X)\times\mathrm{Im}(Y)$, la fibra de un par es
$$(X,Y)^{-1}(\{(x_i,y_j)\})=\{\omega:X(\omega)=x_i\text{ y }Y(\omega)=y_j\}=A_i\cap B_j,$$
pues $X(\omega)=x_i\iff\omega\in A_i$ y $Y(\omega)=y_j\iff\omega\in B_j$. Por tanto la partición es
$$\{A_i\cap B_j:\ 1\le i\le r,\ 1\le j\le q,\ A_i\cap B_j\neq\varnothing\}$$
(solo las intersecciones **no vacías**: un par puede no realizarse). Esto es el **refinamiento
común** de $\{A_i\}$ y $\{B_j\}$: la partición más gruesa que es a la vez más fina que ambas. En
efecto, $A_i=\bigsqcup_j(A_i\cap B_j)$, así que es más fina que $\{A_i\}$, y por simetría más fina
que $\{B_j\}$; y es la más gruesa con esa propiedad. Número de bloques $\le r\cdot q$ y $\le n$.

### Inciso 3 — Con $n=13$ y $X$ fija, ¿qué $Y$ son independientes de $X$?

**Definición.** $X\perp Y\iff P(X=x_i,Y=y_j)=P(X=x_i)P(Y=y_j)$ para todo par.

**Paso 1 — a cardinales.** En el espacio empírico ($P_n(A)=|A|/n$, fibras del inciso 2):
$$\frac{|A_i\cap B_j|}{13}=\frac{|A_i|}{13}\cdot\frac{|B_j|}{13}\iff |A_i\cap B_j|=\frac{|A_i|\,|B_j|}{13}.$$

**Paso 2 — usar que 13 es primo.** El lado izquierdo es un **entero** (un cardinal) ⇒ $13\mid|A_i||B_j|$.
Como 13 es primo, por el **lema de Euclides** $13\mid|A_i|$ **o** $13\mid|B_j|$.

**Paso 3 — acotar.** Cada bloque cumple $1\le|A_i|,|B_j|\le13$, y el único múltiplo de 13 en ese
rango es 13. Luego $13\mid|A_i|\iff|A_i|=13\iff X$ constante; y $13\mid|B_j|\iff|B_j|=13\iff Y$ constante.

**Paso 4 — análisis según $X$ fija:**
- **Caso A — $X$ no constante ($r\ge2$):** todos $|A_i|\le12\Rightarrow13\nmid|A_i|$, así que la
  disyunción obliga $|B_j|=13$ para todo $j\Rightarrow q=1\Rightarrow$ **$Y$ constante**. Recíproca:
  una constante es independiente de cualquier variable (si $|B_1|=13$,
  $|A_i\cap B_1|=|A_i|=\frac{|A_i|\cdot13}{13}$ ✓).
- **Caso B — $X$ constante ($r=1$):** $|A_1|=13$ cumple $13\mid|A_1|$ y la condición vale para
  cualquier $Y$ ⇒ **toda $Y$** es independiente de $X$.

> **Síntesis (papel de la primalidad de 13).** La independencia exige $|A_i\cap B_j|=|A_i||B_j|/n$
> **entero** (estructura de "rectángulo"). Construir dos variables **no triviales** independientes
> requeriría factorizar $n=a\cdot b$ con $1<a,b<n$, **imposible si $n$ es primo**. Por eso sobre
> $\Omega_{13}$ **dos variables solo pueden ser independientes si al menos una es constante**.
> Si $n$ fuera compuesto (p. ej. $n=12=a\times b$), sí existirían variables no constantes
> genuinamente independientes (bloques cuyos tamaños factorizan $n$ en una grilla $a\times b$).
>
> **Respuesta directa al inciso 3:** fijada $X$ — si $X$ es no constante, $Y\perp X\iff Y$ constante;
> si $X$ es constante, toda $Y$ es independiente de $X$.

---

## Apéndice — Chuleta de conceptos transversales

- **Etiqueta vs. variable:** la etiqueta (individuo) **identifica**; la variable **describe**.
  Confundirlas es el error recurrente (preguntas 1, 2, 8).
- **Partición / $\sigma(X)$:** toda variable parte $\Omega_n$ en fibras; $\sigma(X)$ = funciones
  constantes por bloque; $\dim\sigma(X)$ = nº de bloques. Más bloques = más fina.
- **Esperanza condicional = proyección ortogonal** = promedio por bloque = mejor predictor en media cuadrática.
- **Teorema de la Proyección** (existencia/unicidad por completitud + igualdad del paralelogramo)
  ⇒ residuo ortogonal ⇒ $E[\text{residuo}]=0$, esperanza total, Pitágoras, ANOVA.
- **$\mathbf 1\in\sigma(X)$ siempre:** convierte ortogonalidad del residuo en $E[Z]=0$ (insesgamiento).
- **Ortogonalidad lo es todo:** proyecciones se **suman solo si** los subespacios son ortogonales
  (pregunta 3); no correlación ≡ ortogonalidad en $L^2$ centrado (pregunta 7).
- **$R^2$ (centrado) $\le r_n^2$ (crudo):** iguales sii $\mu=0$ o ajuste perfecto (pregunta 4).
- **Norma del supremo** no viene de producto interno (falla el paralelogramo) ⇒ no sirve para
  definir esperanza condicional (pregunta 6).
- **Primalidad y rectángulos:** $n$ primo ⇒ independencia solo con una variable constante (pregunta 8).
- **Frisch–Waugh–Lovell:** para "encajonar" bien con predictores correlacionados, **ortogonalizar primero**.

---

## Tablas de repaso rápido

### Tabla 1 — Fórmulas clave

| Concepto | Fórmula |
|----------|---------|
| Producto interno | $\langle U,V\rangle=E[UV]=\frac1n\sum_{i=1}^n U(\omega_i)V(\omega_i)$ |
| Norma en $L^2$ | $\|U\|^2=E[U^2]=\frac1n\sum_{i=1}^n U^2(\omega_i)$ |
| Esperanza condicional | $E(Y\mid\sigma(X))(\omega)=\frac{1}{|B|}\sum_{\omega'\in B}Y(\omega')$, $B=$ bloque de $\omega$ |
| Proyección sobre recta | $P_{\sigma(v)}Y=\frac{\langle Y,v\rangle}{\langle v,v\rangle}v$ |
| Teorema de Pitágoras | $\|Y\|^2=\|\hat Y\|^2+\|Z\|^2$ |
| Descomposición ANOVA | $\mathrm{Var}(Y)=\mathrm{Var}(E(Y\mid\sigma(X)))+E(\mathrm{Var}(Y\mid\sigma(X)))$ |
| Coef. determinación (crudo) | $r_n^2=\frac{\|E(Y\mid\sigma(X))\|^2}{\|Y\|^2}=\frac{a}{b}$ |
| Coef. determinación (centrado) | $R^2=\frac{\mathrm{Var}(\hat Y)}{\mathrm{Var}(Y)}=\frac{a-\mu^2}{b-\mu^2}$ |
| $r^2$ por bloques (razón corr.) | $r^2=\frac{SC_\text{inter}}{SC_\text{total}}=\frac{\sum_a n_a(\bar x_a-\bar x)^2}{\sum_i x_i^2-n\bar x^2}$ |
| Proyector ortogonal | $P_M x:=\arg\min_{y\in M}\|x-y\|$, única solución con $x-Px\perp M$ |
| Esperanza total | $E[\hat Y]=E[E(Y\mid\sigma(X))]=E[Y]$ |
| Residuo ortogonal | $Z=Y-E(Y\mid\sigma(X))$, $\langle Z,W\rangle=0\ \forall W\in\sigma(X)$ |
| Igualdad del paralelogramo | $\|x+y\|^2+\|x-y\|^2=2\|x\|^2+2\|y\|^2$ (válida en Hilbert) |
| Ecuaciones normales | $X^\top(y-X\beta)=0\Rightarrow\hat\beta=(X^\top X)^{-1}X^\top y$ |
| Matriz proyectora | $P=X(X^\top X)^{-1}X^\top$, $\hat y=Py$, $e=(I-P)y$ |
| Independencia (empírica) | $|A_i\cap B_j|=|A_i||B_j|/n$ para todo $i,j$ |

### Tabla 2 — Símbolos y notación

| Símbolo | Significado |
|---------|-------------|
| $\Omega_n$ | Conjunto finito de etiquetas (individuos), $|\Omega_n|=n$ |
| $\mathcal P(\Omega_n)$ | $\sigma$-álgebra trivial (todas las partes) |
| $P_n$ | Probabilidad uniforme empírica, $P_n(\{\omega\})=1/n$ |
| $\sigma(X)$ | Subespacio de $L^2$ de funciones constantes por bloque de $X$ |
| $L^2(\Omega_n,\mathcal P(\Omega_n),P_n)$ | Espacio de Hilbert de funciones en $\Omega_n$ |
| $\langle U,V\rangle$ | Producto interno en $L^2$: $E[UV]$ |
| $\|U\|$ | Norma de $U$ en $L^2$: $\sqrt{E[U^2]}$ |
| $E(Y\mid\sigma(X))$ | Esperanza condicional = proyección sobre $\sigma(X)$ |
| $\hat Y$ | Predictor óptimo: $E(Y\mid\sigma(X))=P_{\sigma(X)}Y$ |
| $Z$ | Residuo: $Z=Y-\hat Y$ |
| $P_M$ | Proyector ortogonal sobre el subespacio cerrado $M$ |
| $A_i=X^{-1}(\{x_i\})$ | Fibra $i$-ésima de $X$ (bloque con $X=x_i$) |
| $r^2,\eta^2$ | Razón de correlación: $\mathrm{Var}(\hat Y)/\mathrm{Var}(Y)$ |
| $R^2$ | Coef. de determinación (invariante a traslaciones) |
| $r_n^2$ | Coef. basado en momentos crudos (sobrestima si $E[Y]\neq0$) |
| $\mu=E[Y]$ | Esperanza de $Y$ |
| $\mathrm{Var}(Y)$ | Varianza: $E[(Y-\mu)^2]$ |
| $\mathbf 1$ | Función constante $\mathbf 1(\omega)=1$ |
| $M\perp N$ | Subespacios ortogonales |
| $SC_\text{inter/intra/total}$ | Sumas de cuadrados entre / dentro / total |
| $X\perp Y$ | Variables independientes |

### Tabla 3 — Resultados y teoremas

| Teorema / Resultado | Enunciado |
|---------------------|-----------|
| **Caracterización de la proyección** | $Px$ minimiza $\|x-y\|$ sobre $M\iff(x-Px)\perp M$ |
| **Propiedades de $P$** | Lineal, idempotente ($P^2=P$), autoadjunto, contracción ($\|P\|\le1$) |
| **Esperanza = Proyección** | $E[X\mid\mathcal G]=P_M X$, único $\mathcal G$-medible minimizador en $L^2$ |
| **Insesgamiento** | $E[Y-E(Y\mid\sigma(X))]=0$ (porque $\mathbf 1\in\sigma(X)$) |
| **Ortogonalidad del residuo** | $Y-E(Y\mid\sigma(X))\perp\sigma(X)$ |
| **Esperanza total** | $E[E(Y\mid\sigma(X))]=E[Y]$ |
| **Pitágoras** | $\|Y\|^2=\|\hat Y\|^2+\|Y-\hat Y\|^2$ |
| **ANOVA** | $\mathrm{Var}(Y)=\mathrm{Var}(\hat Y)+E[\mathrm{Var}(Y\mid\sigma(X))]$ |
| **$R^2$ vs $r_n^2$** | $R^2\le r_n^2$, iguales sii $E[Y]=0$ o $Y=\hat Y$ c.s. |
| **Suma de proyecciones** | $P_{M\oplus N}=P_M+P_N$ solo si $M\perp N$ |
| **Frisch–Waugh–Lovell** | Regresión parcial: ortogonalizar si los regresores están correlacionados |
| **Coherencia jerárquica** | $\sigma(D)\subseteq\sigma(L)\subseteq\sigma(M)$ (Distrito ⊆ Local ⊆ Mesa) |
| **Falacia ecológica** | Atribuir efecto agregado al nivel individual falla (Robinson, 1950) |
| **Paradoja de Simpson** | La correlación puede invertirse entre nivel agregado e individual |
| **Primalidad e independencia** | Con $n=13$ primo, dos variables no triviales no pueden ser independientes |
| **Igualdad del paralelogramo** | $\|x+y\|^2+\|x-y\|^2=2(\|x\|^2+\|y\|^2)$ caracteriza espacios de Hilbert |
| **Estimador MCO** | $\hat\beta=(X^\top X)^{-1}X^\top y$ es BLUE bajo Gauss–Markov |
| **Ecuaciones normales** | $X^\top e=0\Rightarrow$ residuos ortogonales a regresores |
