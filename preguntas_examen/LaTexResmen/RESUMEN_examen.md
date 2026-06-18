# Resumen para el examen — Estadística geométrica ($L^2$, proyección y esperanza condicional)

> Documento de estudio consolidado a partir de las 8 preguntas de examen.
> Todas comparten **un mismo marco teórico** (enfoque "San Martín"): una base de
> datos es un **espacio empírico de probabilidad**, las variables son **funciones**
> sobre las etiquetas, y la **esperanza condicional es una proyección ortogonal**
> en el espacio de Hilbert $L^2$. Domina este marco y las 8 preguntas se resuelven
> casi mecánicamente.

---

## 0. Marco teórico común (memorízalo: es la llave de todo)

### 0.1 Espacio empírico de probabilidad

Una base de datos es la terna $(\Omega_n, \mathcal{P}(\Omega_n), P_n)$ donde:

- $\Omega_n = \{\omega_1,\dots,\omega_n\}$ es el **conjunto de etiquetas (labels)** = los
  **individuos estadísticos** (lo que identifica a cada fila de forma única). NO confundir
  la etiqueta con el valor de una variable.
- Probabilidad **uniforme** (regla de Laplace): cada etiqueta pesa $p_\omega = 1/n$, y
  $P_n(A) = |A|/n$ para cualquier evento $A\subseteq\Omega_n$.
- Una **variable** es una función $X:\Omega_n \to \mathrm{Im}(X)$. Cada columna de la tabla
  es una variable.

### 0.2 Partición inducida por una variable

Toda variable $X$ con $\mathrm{Im}(X)=\{x_1,\dots,x_r\}$ induce una **partición** de $\Omega_n$
en sus **fibras / conjuntos de nivel / bloques**:
$$A_i := X^{-1}(\{x_i\}) = \{\omega : X(\omega)=x_i\}.$$
Estos bloques son no vacíos, disjuntos dos a dos, y cubren $\Omega_n$ (unión disjunta).
La notación $\sigma(X)$ designa el **subespacio de funciones que son constantes sobre cada bloque**.

### 0.3 $L^2$ como espacio de Hilbert

$L^2(\Omega_n,\mathcal{P}(\Omega_n),P_n)$ con producto interno y norma:
$$\langle U,V\rangle = E[UV] = \sum_{\omega} U(\omega)V(\omega)\,P_n(\{\omega\})
= \tfrac1n\sum_{i=1}^n U(\omega_i)V(\omega_i),\qquad \|U\|^2=\langle U,U\rangle=E[U^2].$$
Como $\Omega_n$ es finito, toda función real está en $L^2$ y el espacio es completo (Hilbert).

### 0.4 Esperanza condicional = proyección ortogonal

$$\boxed{\;E(Y\mid\sigma(X)) = P_{\sigma(X)}\,Y = \text{proyección ortogonal de } Y \text{ sobre } \sigma(X)\;}$$

Es el único elemento $\hat Y\in\sigma(X)$ que minimiza $\|Y-\hat Y\|^2$ (mejor predicción en
media cuadrática). En cada individuo vale **el promedio de $Y$ sobre su bloque**:
$$E(Y\mid\sigma(X))(\omega)=\frac{1}{|B|}\sum_{\omega'\in B} Y(\omega'),\quad B=\text{bloque de }\omega.$$

**Teorema de la Proyección (hecho central, se usa en TODAS las preguntas):** el residuo
$Z=Y-\hat Y$ es ortogonal a todo el subespacio: $\langle Z,W\rangle=0$ para todo $W\in\sigma(X)$.

**Dato clave:** la función constante $\mathbf 1$ siempre está en $\sigma(X)$ (es constante en
cada bloque). De ahí salen casi todas las identidades de esperanza.

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
  $P_{\mathcal M\oplus\mathcal N}=P_{\mathcal M}+P_{\mathcal N}$. Si no son ortogonales, esto **falla**
  (hay que resolver ecuaciones normales).

---

## 1. Base de datos (impacto de reformas fiscales)

**Planteo:** tabla País×Año (20 individuos, 9 variables). Espacio $\mathbb R^{20}$, pesos $1/20$.

| # | Pregunta | Respuesta clave |
|---|----------|-----------------|
| 1 | Conjunto de etiquetas | Los **pares (País, Año)**: $I=\{(\mathrm{FRA},2015),\dots\}$, $|I|=4\times5=20$. Cada label pesa $1/20$. (El enunciado nombra 5 países pero la tabla solo trae 4.) |
| 2 | Partición inducida por FR | FR es binaria $\Rightarrow$ **2 clases**: $A_0=\{FR=0\}$ (2015–2016, sin reforma, $|A_0|=8$) y $A_1=\{FR=1\}$ (2017–2019, con reforma, $|A_1|=12$). |
| 3 | $E(\text{Inv}\mid\text{Año})$ | Año $\Rightarrow$ 5 clases (4 individuos c/u). En cada individuo vale el **promedio del año**: 2015→18.175, 2016→17.300, 2017→18.825, 2018→21.000, 2019→18.325. Media global $\bar x=374.5/20=18.725$. |
| 4 | $r^2$ e interpretación geométrica | Razón de correlación $\eta^2$. |

**Cálculo del $r^2$ (pregunta 4):**
$$SC_\text{inter}=\sum_a n_a(\bar x_a-\bar x)^2 = 30.715,\quad
SC_\text{total}=\sum_i x_i^2-20\bar x^2 = 53.8975,\quad SC_\text{intra}=23.1825.$$
$$r^2=\frac{SC_\text{inter}}{SC_\text{total}}=\frac{30.715}{53.8975}\approx 0.570 \;(\approx 57\%).$$

**Interpretación geométrica:** $E(\text{Inv}\mid\text{Año})$ es la proyección ortogonal de Inv
sobre el subespacio $F$ de funciones constantes por año; el residuo es ortogonal a $F$. Con
vectores centrados, $\|Inv-\bar x\|^2=\|\hat Y-\bar x\|^2+\|Inv-\hat Y\|^2$ (Pitágoras) y
$r^2=\cos^2\theta$, donde $\theta$ es el ángulo entre Inv (centrada) y $F$. Aquí
$\cos\theta=\sqrt{0.57}\approx0.755\Rightarrow\theta\approx41^\circ$.

---

## 2. Coherencia en una base de datos (datos electorales / falacia ecológica)

**Planteo:** votos identificados por folio único; cada folio tiene opción de voto y pertenece a
mesa → local → distrito.

1. **Etiquetas = los folios** (un voto = un individuo). La **opción de voto NO es etiqueta**, es una
   variable $V:\Omega_N\to\{\text{cand.1, cand.2, nulo, blanco}\}$. Error conceptual de fondo:
   confundir la etiqueta con el valor de una variable. Los folios deben ser únicos, exhaustivos y
   mutuamente excluyentes.
2. **Mesa, Local, Distrito SÍ son variables** (funciones bien definidas sobre los folios). Sus
   particiones son **jerárquicas/encajonadas**: Mesa refina a Local que refina a Distrito:
   $$\sigma(D)\subseteq\sigma(L)\subseteq\sigma(M)\subseteq L^2(\Omega_N).$$
   (Mesa = más fina; Distrito = más gruesa. Existe función de paso mesa→local→distrito, no al revés.)
3. **Condiciones para incorporar la tasa de cesantía comunal $C$:**
   - **Coherencia jerárquica:** debe existir asignación unívoca folio→comuna (cada local en una sola
     comuna), definiendo $K:\Omega_N\to\{\text{comunas}\}$ compatible con la jerarquía.
   - **Medibilidad:** $C$ debe ser constante por comuna, i.e. $C\in\sigma(K)$.

**El punto de fondo (falacia ecológica):** como $C\in\sigma(K)$ es constante dentro de cada comuna,
solo aporta variabilidad **entre** comunas, no dentro. Toda regresión factoriza por
$E(V\mid C)=E(V\mid\sigma(K))$ y no distingue dos votantes de la misma comuna:
$$\mathrm{Var}(V)=\underbrace{\mathrm{Var}(E(V\mid\sigma(K)))}_{\text{entre comunas (lo que ve }C)}
+\underbrace{E(\mathrm{Var}(V\mid\sigma(K)))}_{\text{dentro: invisible para }C}.$$
Atribuir el efecto agregado al individuo es la **falacia ecológica** (Robinson, 1950); puede incluso
invertirse (**paradoja de Simpson**). Para explicar el nivel individual hace falta una variable
definida sobre el folio con partición más fina que la comunal.

---

## 3. Modelos lineales encajonados (nested) — crítica vía ortogonalidad

**Modelos:** M1: $E(Y\mid X_1)=\alpha_0+\alpha_1 X_1$; M2: $E(Y\mid X_1,X_2)=\alpha_0+\alpha_1X_1+\alpha_2X_2$.
Práctica habitual: "ajusto M1, agrego $X_2$ y veo si mejora". Los ejemplos critican cuándo es válida.

**Marco:** $\Omega_3$, prob. uniforme $1/3$, vectores en $\mathbb R^3$, $\langle U,V\rangle=\tfrac13\sum U_iV_i$.

**Ejemplo 1 — regresores ORTOGONALES (se cumple la aditividad).**
$\mathbf 1=(1,1,1)$, $X=(-1,0,1)$. Como $\langle\mathbf 1,X\rangle=\tfrac13(-1+0+1)=0$, son ortogonales:
$$E(Y\mid\sigma(\mathbf 1,X))=\bar y\,\mathbf 1+\tfrac{y_3-y_1}{2}X
= E(Y\mid\sigma(\mathbf 1))+E(Y\mid\sigma(X)).\;\checkmark$$
**Causa de la igualdad:** ortogonalidad de los regresores ⇒ proyectar sobre el plano = sumar las
proyecciones sobre cada eje (Pitágoras / descomposición ortogonal). No hay "magia".

**Ejemplo 2 — regresores NO ortogonales (falla la aditividad).**
$\mathbf 1=(1,1,1)$, $Z=(1,0,1)$. Ahora $\langle\mathbf 1,Z\rangle=\tfrac23\neq0$. Hay que resolver
**ecuaciones normales** ($a\langle\mathbf 1,\mathbf 1\rangle+b\langle Z,\mathbf 1\rangle=\langle Y,\mathbf 1\rangle$, etc.):
$$a=y_2,\quad b=\tfrac{y_1+y_3}{2}-y_2,\quad
E(Y\mid\sigma(\mathbf 1,Z))=\Big(\tfrac{y_1+y_3}{2},\,y_2,\,\tfrac{y_1+y_3}{2}\Big).$$
Comparando en $\omega_2$: $\big[E(Y\mid\sigma(\mathbf 1,Z))\big]=y_2$ vs.
$\big[E(Y\mid\sigma(\mathbf 1))+E(Y\mid\sigma(Z))\big]=\tfrac{y_1+y_2+y_3}{3}$. Coinciden solo si
$2y_2=y_1+y_3$ (caso particular). En general:
$$E(Y\mid\sigma(\mathbf 1,Z))\neq E(Y\mid\sigma(\mathbf 1))+E(Y\mid\sigma(Z)).$$
Al no ser ortogonales, $\sigma(\mathbf 1)$ y $\sigma(Z)$ se **solapan**; sumar proyecciones cuenta dos
veces la parte común.

**Respuesta (¿cuándo son genuinamente "encajonados"?):** Hay que distinguir dos nociones:
1. **Inclusión de subespacios:** $\sigma(\mathbf 1)\subseteq\sigma(\mathbf 1,X_1)\subseteq\sigma(\mathbf 1,X_1,X_2)$
   **siempre** se cumple. Solo valida comparar ajustes (M2 nunca ajusta peor).
2. **Descomposición aditiva de la contribución** (lo que la práctica realmente usa): se cumple
   **si y solo si el nuevo predictor es ortogonal a los anteriores**, $X_2\perp\sigma(\mathbf 1,X_1)$.

> **Criterio:** M1 y M2 son "encajonados" en sentido fuerte (el coef. de $X_1$ no cambia al agregar
> $X_2$ y la mejora se atribuye limpiamente a $X_2$) **solo cuando $X_2\perp\sigma(\mathbf 1,X_1)$**.
> Si los predictores están correlacionados, agregar $X_2$ altera el coef. de $X_1$. La forma correcta
> de "encajonar" es **ortogonalizar primero** (idea del **teorema de Frisch–Waugh–Lovell**): comparar
> $Y$ contra la parte de $X_2$ residual respecto de $\sigma(\mathbf 1,X_1)$. La receta ingenua
> "ajusto, agrego, comparo" solo es interpretable bajo ortogonalidad.

---

## 4. Teorema de Pitágoras ($R^2$ vs $r_n^2$)

**Setup:** $Y=\hat Y+Z$ con $\hat Y=E(Y\mid\sigma(X))$, $Z=Y-\hat Y$, y $\|Y\|^2=\|\hat Y\|^2+\|Z\|^2$.

**Parte 1 — $E(Y-E(Y\mid\sigma(X)))=0$:** como $\mathbf 1\in\sigma(X)$ y $Z\perp\sigma(X)$,
$\langle Z,\mathbf 1\rangle=E[Z]=0$. Consecuencia: $E[\hat Y]=E[Y]=:\mu$ (insesgamiento).

**Parte 2 — equivalencia (1)⟺(2):** usando $\mathrm{Var}(U)=E[U^2]-[E(U)]^2$ y que $E[\hat Y]=E[Y]=\mu$,
$E[Z]=0$. Restando $\mu^2$ a la identidad de Pitágoras $E[Y^2]=E[\hat Y^2]+E[Z^2]$:
$$\mathrm{Var}(Y)=\mathrm{Var}(E(Y\mid\sigma(X)))+\mathrm{Var}(Y-E(Y\mid\sigma(X))).$$
La operación "restar $\mu^2$" es reversible ⇒ (1) y (2) son equivalentes. La equivalencia descansa en
la Parte 1 (el ajuste por la media $-\mu^2$ es común).

**Parte 3 — comparar $R^2$ y $r_n^2$:** con $a=\|\hat Y\|^2=E[\hat Y^2]$, $b=\|Y\|^2=E[Y^2]$:
$$r_n^2=\frac{\|E(Y\mid\sigma(X))\|^2}{\|Y\|^2}=\frac{a}{b}
\qquad\text{(segundos momentos crudos)},$$
$$R^2=\frac{\mathrm{Var}(\hat Y)}{\mathrm{Var}(Y)}=\frac{a-\mu^2}{b-\mu^2}
\qquad\text{(varianzas, centrado en la media)}.$$
Ordenamiento: $b=a+E[Z^2]\Rightarrow a\le b$; por Jensen $a=E[\hat Y^2]\ge\mu^2$. Luego $\mu^2\le a\le b$.
Producto cruzado ⇒
$$\boxed{R^2\le r_n^2.}$$
**Igualdad** sii $\mu^2(b-a)=0$, i.e. $\mu=E[Y]=0$ (variable centrada) o $a=b$ ($Y=\hat Y$ c.s., ajuste perfecto).

**Interpretación:** $R^2$ es el **genuino coeficiente de determinación** (proporción de varianza
explicada, invariante a traslaciones); $r_n^2$ usa momentos crudos, incorpora la media y **sobrestima**
el ajuste. La diferencia crece con $|\mu|$.

**Ejemplo concreto:** $\Omega=\{1,2,3,4\}$ uniforme, bloques $\{1,2\},\{3,4\}$, $Y=(0,2,2,4)$.
$\hat Y=(1,1,3,3)$, $\mu=2$, $b=6$, $a=5$, $\mathrm{Var}(Y)=2$, $\mathrm{Var}(\hat Y)=1$, $\mathrm{Var}(Z)=1$.
$$r_n^2=\tfrac{5}{6}\approx0.833,\qquad R^2=\tfrac{1}{2}=0.5.\quad (R^2<r_n^2\;\checkmark)$$

---

## 5. Esperanza condicional (construcción de la proyección + ley total)

**Parte 1 — construir $P:H\to M$ y probar que es proyección ortogonal** ($M\subseteq H$ cerrado):

- **Definición:** $Px:=y_0$, el único minimizador de $\|x-y\|$ sobre $M$ (bien definido por unicidad).
- **Lema clave (caracterización):** $y_0$ minimiza $\iff (x-y_0)\perp M$.
  - ($\Leftarrow$) Si $x-y_0\perp M$, Pitágoras da $\|x-y\|^2=\|x-y_0\|^2+\|y_0-y\|^2\ge\|x-y_0\|^2$.
  - ($\Rightarrow$) Si minimiza, $\varphi(t)=\|x-(y_0+tm)\|^2$ tiene mínimo en $t=0\Rightarrow\varphi'(0)=-2\langle x-y_0,m\rangle=0$.
- **Lineal:** $(ax_1+bx_2)-(aPx_1+bPx_2)\perp M$ ⇒ por unicidad $P(ax_1+bx_2)=aPx_1+bPx_2$.
- **Idempotente:** $P|_M=\mathrm{id}\Rightarrow P^2=P$.
- **Autoadjunto:** $\langle Px,z\rangle=\langle Px,Pz\rangle=\langle x,Pz\rangle$.
- **Contracción y descomposición:** $\|Px\|\le\|x\|$ ($\|P\|\le1$) y $H=M\oplus M^\perp$.

**Parte 2 — por qué la proyección ES la esperanza condicional:** en $H=L^2(\Omega,\mathcal F,P)$ con
$\langle X,Y\rangle=E[XY]$, el subespacio $M=L^2(\Omega,\mathcal G,P)$ ($\mathcal G$-medibles) es cerrado.
La definición de $E[X\mid\mathcal G]$ (la variable $\mathcal G$-medible $Z$ con $E[(X-Z)W]=0\;\forall W\in M$)
es **exactamente** la condición de ortogonalidad del lema ⇒ $E[X\mid\mathcal G]=PX$. Interpretación:
mejor predicción de $X$ en media cuadrática usando solo la información de $\mathcal G$.

**Parte 3 — demostrar $E[E(Y\mid\sigma(X))]=E(Y)$:** la función constante $\mathbf 1\in M=\sigma(X)$, y
$E[Z]=\langle Z,\mathbf 1\rangle$. Como $Y-\hat Y\perp M$ y $\mathbf 1\in M$:
$$E[Y]-E[\hat Y]=\langle Y-\hat Y,\mathbf 1\rangle=0\;\Rightarrow\;E[\hat Y]=E[Y].$$
Verificación directa: $\hat Y|_{B_j}=E[Y\mid B_j]$ y $E[\hat Y]=\sum_j P(B_j)E[Y\mid B_j]=E[Y]$ (**ley de la esperanza total**).

---

## 6. Igualdad del paralelogramo

**Enunciado:** en un espacio con producto interno, $\|x+y\|^2+\|x-y\|^2=2\|x\|^2+2\|y\|^2$.

**1. Demostración:** expandir con $\|v\|^2=\langle v,v\rangle$ (única vez que se usa que la norma viene de
un producto interno) usando bilinealidad y simetría:
$$\|x+y\|^2=\|x\|^2+2\langle x,y\rangle+\|y\|^2,\qquad \|x-y\|^2=\|x\|^2-2\langle x,y\rangle+\|y\|^2.$$
Al sumar, los términos cruzados $\pm2\langle x,y\rangle$ **se cancelan** ⇒ resultado. (Caso complejo:
simetría hermitiana, cruzados $=2\,\mathrm{Re}\langle x,y\rangle$, cancela igual.)

**2. Relevancia para la esperanza condicional:** la esperanza condicional es una **proyección ortogonal
en $L^2$**, y su existencia/unicidad es el **Teorema de la Proyección**, cuya prueba usa la igualdad del
paralelogramo de forma esencial: dada una sucesión minimizante $(Y_n)$ con $\|X-Y_n\|\to d$, aplicando la
igualdad a $a=X-Y_n$, $b=X-Y_m$:
$$\|Y_m-Y_n\|^2=2\|X-Y_n\|^2+2\|X-Y_m\|^2-4\Big\|X-\tfrac{Y_n+Y_m}{2}\Big\|^2\le 2\|X-Y_n\|^2+2\|X-Y_m\|^2-4d^2\to0,$$
luego $(Y_n)$ es de **Cauchy**; por completitud de $L^2$ converge a $Y^*\in M$ (cerrado), el minimizador.
Sin la igualdad del paralelogramo (norma sin producto interno) el argumento se rompe ⇒ por eso la
esperanza condicional se define en $L^2$ y no con la norma del supremo.

**3. La norma del supremo en $\mathbb R^n$ ($n\ge2$) NO la satisface:** $x=e_1$, $y=e_2$.
$\|x\|_\infty=\|y\|_\infty=1$, $\|x+y\|_\infty=\|x-y\|_\infty=1$. LI $=1^2+1^2=2$; LD $=2\cdot1+2\cdot1=4$.
$2\neq4$ ⇒ falla ⇒ no proviene de un producto interno.

**4. La norma del supremo en $C[0,1]$ NO la satisface:** $F(x)=1$, $G(x)=x$.
$\|F\|_\infty=\|G\|_\infty=1$, $\|F+G\|_\infty=\sup|1+x|=2$, $\|F-G\|_\infty=\sup|1-x|=1$.
LI $=2^2+1^2=5$; LD $=4$. $5\neq4$ ⇒ falla.

---

## 7. El término de error — dos miradas

**Modelo:** $y=X\beta+\varepsilon$, con $X\in\mathbb R^{n\times p}$, $\mathrm{rango}(X)=p$ (⇒ $X^\top X$ invertible).
Ambas miradas llegan al **mismo estimador** $\hat\beta=(X^\top X)^{-1}X^\top y$ por caminos opuestos.

### Mirada de Gauss (probabilística)
- $\varepsilon$ es **objeto primitivo**: variable aleatoria con ley **supuesta**.
- Supuestos Gauss–Markov: $E[\varepsilon]=0$ (no sesga), $\mathrm{Var}(\varepsilon)=\sigma^2 I_n$
  (**homocedasticidad** + **no correlación**). Para inferir: $\varepsilon\sim\mathcal N(0,\sigma^2 I_n)$.
- **Argumento histórico:** la media aritmética como estimador más plausible ⇒ la distribución que la hace
  máximo-verosímil es la **normal** ⇒ maximizar verosimilitud normal = **minimizar $\sum\varepsilon_i^2$**.
  Después (Theoria Combinationis): **Gauss–Markov** ⇒ sin normalidad, $\hat\beta$ es **BLUE** (mínima
  varianza entre lineales insesgados).
- Estimación: minimizar $S(\beta)=\|y-X\beta\|^2$, derivar ⇒ $X^\top X\hat\beta=X^\top y$.
- **Estilo:** paramétrico, inferencial, muestral (insesgadez, $\mathrm{Var}(\hat\beta)=\sigma^2(X^\top X)^{-1}$,
  tests $t$ y $F$, IC, Cramér–Rao). Mentalidad: "de la probabilidad hacia las propiedades".

### Mirada del Teorema de la Proyección (geométrica)
- Punto de partida: $\mathbb R^n$ con producto interno, $\mathcal C(X)=\{X\beta\}$.
- $y=\hat y+e$ con $\hat y\in\mathcal C(X)$, $e\perp\mathcal C(X)$. El error **se define** como complemento
  ortogonal (no se supone su ley, se **construye**).
- **Ecuaciones normales = ortogonalidad** (sin derivar): $X^\top e=0\iff X^\top X\hat\beta=X^\top y$.
  Matriz de proyección $P=X(X^\top X)^{-1}X^\top$ ($P^2=P$, $P^\top=P$), $\hat y=Py$, $e=(I-P)y$.
- **Pitágoras = ANOVA sin distribución:** $\|y\|^2=\|\hat y\|^2+\|e\|^2$ ⇒ $SCT=SCR+SCE$.
- **Versión profunda en $L^2$:** ortogonalidad ≡ no correlación; $E[Y\mid X]$ es la proyección de $Y$;
  el error $\varepsilon=Y-E[Y\mid X]$ es **automáticamente** no correlacionado con los regresores
  (lo que en Gauss era hipótesis, aquí es teorema gratuito).
- **Estilo:** geométrico, estructural, libre de distribución. Mentalidad: "de la estructura del espacio
  hacia las ecuaciones; la probabilidad se agrega después, solo para inferir".

### Síntesis
| Criterio | Gauss | Proyección |
|---|---|---|
| Naturaleza del error | Primitivo (v.a. con ley) | Derivado (residuo ortogonal) |
| No correlación | **Supuesto** | **Consecuencia** automática |
| Ecuaciones normales | Anular $\partial S/\partial\beta$ | Traducir $X^\top e=0$ |
| Por qué MCO | Verosimilitud normal / BLUE | Proyección = punto más cercano |
| Distribución | Central desde el inicio | Postergada (solo para inferir) |
| Estilo | Inferencial/paramétrico | Geométrico/estructural |

Gauss responde **por qué es óptimo** (aporte inferencial); la proyección responde **por qué tiene esa
forma** (aporte estructural). Distinción fina: **errores verdaderos** $\varepsilon$ ≠ **residuos** $e$;
$e=(I-P)\varepsilon$ es una estimación de $\varepsilon$. La síntesis moderna combina ambas: geometría
para la forma/propiedades, probabilidad para la inferencia.

---

## 8. Independencia de funciones ($n=13$ primo)

**Setup:** $\Omega_n$ empírico uniforme, $X:\Omega_{13}\to\{x_1,\dots,x_r\}$, $Y:\Omega_{13}\to\{y_1,\dots,y_q\}$.

**1. ¿$r$ arbitrario?** No: $\boxed{1\le r\le n}$. Como los bloques $A_i$ son no vacíos y particionan,
$n=\sum_{i=1}^r|A_i|\ge\sum_{i=1}^r 1=r$. Casos extremos: $r=1\iff X$ constante; $r=n\iff X$ inyectiva.
Con $n=13$: $r\in\{1,\dots,13\}$.

**2. Partición inducida por $(X,Y)$:** la fibra de un par es
$$(X,Y)^{-1}(\{(x_i,y_j)\})=A_i\cap B_j,$$
así que la partición es $\{A_i\cap B_j: A_i\cap B_j\neq\varnothing\}$ (solo intersecciones no vacías). Es el
**refinamiento común** (la partición más gruesa que es a la vez más fina que $\{A_i\}$ y que $\{B_j\}$).
Número de bloques $\le r\cdot q$ y $\le n$.

**3. Con $n=13$ y $X$ fija, ¿qué $Y$ son independientes?** Independencia:
$P(X=x_i,Y=y_j)=P(X=x_i)P(Y=y_j)$. En el espacio empírico ($P_n(A)=|A|/n$):
$$\frac{|A_i\cap B_j|}{13}=\frac{|A_i|}{13}\cdot\frac{|B_j|}{13}\iff |A_i\cap B_j|=\frac{|A_i||B_j|}{13}.$$
El lado izquierdo es entero ⇒ $13\mid |A_i||B_j|$. Como **13 es primo** (lema de Euclides):
$13\mid|A_i|$ o $13\mid|B_j|$. Pero $1\le|A_i|,|B_j|\le13$, y el único múltiplo de 13 en ese rango es 13:
$$13\mid|A_i|\iff |A_i|=13\iff X \text{ constante};\qquad 13\mid|B_j|\iff Y\text{ constante}.$$

- **Caso A ($X$ no constante, $r\ge2$):** todos $|A_i|\le12$ ⇒ $13\nmid|A_i|$ ⇒ se necesita $|B_j|=13$ para
  todo $j$ ⇒ $q=1$ ⇒ **$Y$ constante**. Recíproca: una constante es independiente de cualquier variable.
- **Caso B ($X$ constante):** **toda $Y$** es independiente de $X$.

> **Síntesis (papel de la primalidad de 13):** la independencia exige $|A_i\cap B_j|=|A_i||B_j|/n$ entero
> (estructura de "rectángulo"); construir dos variables no triviales independientes requeriría factorizar
> $n=a\cdot b$ con $1<a,b<n$, imposible si $n$ es **primo**. Por eso sobre $\Omega_{13}$ **dos variables solo
> pueden ser independientes si al menos una es constante**. Si $n$ fuera compuesto (p. ej. $n=12=a\times b$)
> sí existirían variables no constantes genuinamente independientes.

---

## Apéndice — Chuleta de conceptos transversales

- **Etiqueta vs. variable:** la etiqueta (individuo) identifica; la variable describe. Confundirlas es el
  error recurrente (preguntas 1, 2, 8).
- **Partición / $\sigma(X)$:** toda variable parte $\Omega_n$ en fibras; $\sigma(X)$ = funciones constantes
  por bloque. Refinamiento: más bloques = más fina.
- **Esperanza condicional = proyección ortogonal** = promedio por bloque = mejor predictor en media cuadrática.
- **Teorema de la Proyección** (existencia/unicidad por completitud + igualdad del paralelogramo) ⇒ residuo
  ortogonal ⇒ $E[\text{residuo}]=0$, esperanza total, Pitágoras, ANOVA.
- **Ortogonalidad lo es todo:** proyecciones se suman solo si los subespacios son ortogonales (pregunta 3);
  no correlación ≡ ortogonalidad en $L^2$ centrado (pregunta 7).
- **$R^2$ (centrado) $\le r_n^2$ (crudo)**; iguales si $\mu=0$ o ajuste perfecto (pregunta 4).
- **Norma del supremo** no viene de producto interno (falla paralelogramo) ⇒ no sirve para definir
  esperanza condicional (pregunta 6).
