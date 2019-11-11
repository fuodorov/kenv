
# Огибающая пучка с учетом влияния пространственного заряда

<a href=mailto:fuodorov1998@gmail.com>V. Fedorov</a>, <a href=mailto:nikdanila@bk.ru>D. Nikiforov</a>, <a href=http://www.inp.nsk.su/~petrenko/>A. Petrenko</a>, (Novosibirsk, 2019)


Сначала разберем теорию, затем расчитаем огибающую электронного пучка в ЛИУ-5 с помощью Python и Astra, приведем сравнение результатов

### Уравнения Максвелла

Объемная плотность тока пучка $\jmath = \rho\upsilon$, где $\rho$ - объемная плотность заряда,$\upsilon$ - скорость пучка.
Запишем дифференциальные уравнения Максвелла:
$$
\nabla \vec{D} = 4\pi\rho,
$$
$$
\nabla\times \vec{H} = \frac{4\pi\vec{\jmath}}{c},
$$
где $\vec{D}$ - индукция электрического поля, $\vec{H}$ - напряженность магнитного поля, $с$ - скорость света. Используем теорему Стокса об интегрировании дифференциальных форм, чтобы получить уравнения Максвелла в интегральной форме:
$$
\oint\limits_{\eth V} \vec{D}\vec{dS} = 4\pi\int\limits_V\rho{dV},
$$
$$
\oint\limits_{\eth S} \vec{H}\vec{dl} = \frac{4\pi}{c}\int\limits_S\vec{\jmath}\vec{dS}.
$$

Найдем $D_r$ для цилиндрического пучка радиуса $a$ с постоянной плотностью $\rho_0$:
$$
D_r = \frac{4\pi}{r}\int\limits^r_0\rho(\xi)\xi d\xi = 
\begin{equation*}
 \begin{cases}
   \displaystyle 2\pi\rho_0 r, r < a, 
   \\
   \displaystyle \frac{2\pi\rho_0 a^2}{r}, r > a.
 \end{cases}
\end{equation*}
$$

Учитывая, что в вакууме 
$D = E$, $H = B$, $E$ - напряженность электрического поля, 
$B$ - индукция магнитного поля, и в плоском пространстве в декартовой системе координат 
$H_\alpha = \beta D_r$, 
где $\displaystyle\beta = \frac{\upsilon}{c}$, 
радиальная компонента силы $F_r$ из силы Лоренца 
$\vec{F} = e\vec{E} + \displaystyle\frac{e}{c}\vec{\upsilon}\times\vec{B}$ :
$$
F_r = eE_r - \frac{e\upsilon_z B_\alpha}{c} = eE_r(1-\displaystyle \frac{\upsilon^2}{c^2}) = \displaystyle \frac{eE_r}{\gamma^2}.
$$
Полезно выразить поле через ток
$
I = \rho_0\upsilon\pi a^2
$, тогда: 
$$
E_r = 
\begin{equation*}
 \begin{cases}
   \displaystyle \frac{2 I r}{a^2\upsilon}, r < a, 
   \\
   \displaystyle \frac{2 I}{r\upsilon}, r > a.
 \end{cases}
\end{equation*}
$$

### Уравнения движения

Второй закон Ньютона $\dot p_r = F_r$, используем параксиальное приближение и считаем $\gamma = const$:
$$
\gamma m \ddot r = \displaystyle \frac{eE_r}{\gamma^2} = \displaystyle \frac{2 I e}{\gamma^2 a^2 \upsilon} r, 
$$
получилось линейное уравнение, но так как все частицы движутся нужно учесть, что $a = a(t)$.
Решение линейного уравнения можно представить как линейное преобразование фазовой плоскости. Так как отрезок на фазовой плоскости при невырожденном линейном преобразовании переходит в отрезок, то его можно охарактеризовать одной точкой. Следовательно, выберем крайнюю точку $r = a$, которая будет характеризовать крайнюю траекторию:
$$
\gamma m \ddot r = \displaystyle \frac{2 I e}{\gamma^2 a \upsilon}. 
$$
Перейдем к дифференцированию по $z$, учтем, что $\displaystyle dt = \frac{dz}{v}$, тогда:
$$
a'' = \displaystyle \frac{e}{a}\frac{2I}{m\gamma^3\upsilon^3}.
$$
Введем характерный альфвеновский ток $I_a = \displaystyle \frac{mc^3}{e} \approx$ 17 kA, следовательно:
$$
a'' = \displaystyle \frac{2I}{I_a (\beta\gamma)^3} \frac{1}{a}.
$$
Учтем внешнюю фокусировку, предполагая суперпозицию полей (верно не всегда, например, в нелинейных средах это не выполняется), получим:
$$
a'' + k(z)a - \displaystyle \frac{2I}{I_a (\beta\gamma)^3} \frac{1}{a} = 0,
$$
что напоминает уравнение огибающей: 
$$ w'' + kw - \displaystyle \frac{1}{w^3} = 0 ,$$
где $ w =\displaystyle  \sqrt \beta .$

### Уравнения огибающей для эллиптического пучка с распределением Капчинского-Владимирского с внешней фокусировкой линейными полями

Распределение Капчинского-Владимирского:
$$
f = A\delta(1 - \displaystyle\frac{\beta_x x'^2 + 2\alpha_x x x' + \gamma_x x^2}{\epsilon_x} - \displaystyle\frac{\beta_y y'^2 + 2\alpha_y y y' + \gamma_y y^2}{\epsilon_y} ),
$$
где $А$ - инвариант Куранта-Снайдера. Полуоси эллипса:
$$
a = \sqrt{\epsilon_x \beta_x}, b = \sqrt{\epsilon_y \beta_y}. 
$$
Поле получается линейно внутри заряженного эллиптического цилиндра:
$$
E_x = \displaystyle \frac{4I}{\upsilon}\frac{x}{a(a+b)},
$$
$$
E_y = \displaystyle \frac{4I}{\upsilon}\frac{y}{b(a+b)}.
$$
Проверим, что $\nabla \vec{E} = 4\pi\rho:$
$$
\displaystyle I = \rho \upsilon \pi ab,
$$
$$
\nabla \vec{E} = \displaystyle \frac{4I(a+b)}{\pi(a+b)ab} = \displaystyle \frac{4I}{\pi ab} = 4\pi\rho.
$$
Так как поля линейные, они добавятся к полям фокусирующей линзы. Подставим в уравнение огибающей $\displaystyle  a = \sqrt \epsilon_x w_x, b = \sqrt \epsilon_y w_y:$
$$
a'' + k_{xt} a - \frac{\epsilon_x^2}{a^3} = 0,
$$
где $k_{xt} = k_x + k_{xsc}$ - полная жесткость, $k_x$ - жесткость линзы, а $\displaystyle k_{xsc} = \frac{4I}{I_a (\beta\gamma)^3}\frac{1}{a(a+b)}.$
В итоге получаем систему уравнений, связанных через пространственный заряд:
$$
\begin{equation*}
 \begin{cases}
   \displaystyle a'' + k_xa - \frac{4I}{I_a (\beta\gamma)^3}\frac{1}{(a+b)} - \frac{\epsilon_x^2}{a^3} = 0 ,
   \\
   \displaystyle b'' + k_yb - \frac{4I}{I_a (\beta\gamma)^3}\frac{1}{(a+b)} - \frac{\epsilon_y^2}{b^3} = 0.
 \end{cases}
\end{equation*}
$$


### Количественный критерий применимости приближения ламинарности течения

Данная система уравнений позволяет учесть 2 эффекта, мешающих сфокусировать пучок в точку - конечность эммитанса и пространственный заряд. Работают члены одинаково, поэтому можно сравнить эти величины.
Когда ток $I$ малый - слабое отталкивание, если ток $I$ большой - сильное отталкивание, следовательно, эммитанс можно откинуть и считать течение ламинарным. Очевидно, количественный критерий применимости ламинарности течения выглядит так: 
$$
\displaystyle \sqrt{\frac{2I}{I_a(\beta\gamma)^3}} \gg \sqrt{\frac{\epsilon}{\beta}}.
$$
Видно, что пространственный заряд влияет на огибающую больше там, где $\beta$-функция больше, а в вблизи фокуса влиянием пространственного заряда можно пренебречь.

### Теорема Буша

В качестве дополнения к выводу уравнения огибающей пучка докажем важную вспомогательную теорему, которая называется теоремой Буша. Она связывает угловую скорость заряженной частицы, движущейся в аксиально-симметричном магнитном поле, с магнитным потоком, охваченным окружностью с центром на оси и проходящим через точку, в которой расположена частица.

Рассмотрим заряд $q$, движущийся в магнитном поле $\vec B = (B_r,0,B_z)$. Приравняем $\theta $ -составляющую силы Лоренца к производной момента импульса по времени, деленной на $r$:
$$
F_\theta = -q(\ddot r B_z - \dot z B_r) = \frac{d}{rdt}(\gamma m r^2 \dot \theta).
$$
Поток, пронизывающий площадь, охваченную окружностью радиуса $r$, центр которой расположен на оси, а сама она проходит через точку, в которой расположен заряд, записывается в виде $\psi = \int\limits^r_0 2\pi r B_z dr$. Когда частица перемещается на $\vec{dl} = (dr,dz)$, скорость изменения потока, охваченного этой окружностью, можно найти из второго уравнения Максвелла $\nabla \vec{B} = 0.$ Таким образом,
$$
\dot\psi = 2\pi r (-B_r \dot z + B_z \dot r).
$$
После интегрирования по времени из приведенных уравнений получаем следующее выражение:
$$
\dot \theta = (-\frac{q}{2\pi \gamma m r^2})(\psi - \psi_0).
$$

### Уравнение параксиального луча

Здесь мы запишем уравнение параксиального луча в виде, соответсвующем системе с аксиальнорй симметрией при принятых ранее допущениях. Чтобы вывести уравнение параксиального луча, приравняем силу радиального ускорения силам электрическим и магнитным со стороны внешних полей. Нужно помнить, что величина $\gamma = \gamma(t).$
$$
\frac{d}{dt}(\gamma m \dot r) - \gamma m r (\dot \theta)^2 = q(E_r + r\dot \theta B_z).
$$
Применим теорему Буша и независимость $B_z$ от $r:$
$$
 -\dot \theta = \frac{q}{2\gamma m}(B_z - \frac{\psi_0}{\pi r}).
$$
Исключим $\dot \theta$ и подставим $\displaystyle \dot \gamma \approx \frac {\beta q E_z }{mc}:$
$$
\ddot r + \frac{\beta q E_z}{\gamma m c} \dot r + \frac{q^2 B^2_z}{4\gamma ^2 m^2} r - \frac{ q^2 \psi^2_0}{4\pi^2\gamma^2 m^2} (\frac{1}{r^3}) - \frac{q E_r}{\gamma m } = 0.
$$

### Уравнение огибающей для круглого и эллиптического пучка

Учитывая, что:
$$
\dot r = \beta c r',
$$
$$
\ddot r = r'' (\dot z)^2 + r'\ddot z \approx r'' \beta^2 c^2 + r' \beta' \beta c^2.
$$
А также, если в области пучка нет никаких зарядов, то, разлагая в ряд Тейлора в окрестности оси и оставляя только первый член, с учетом $$\nabla \vec{E} = 0$$ получаем
$$
E_r = -0.5 r E'_z \approx - 0.5 r \gamma'' mc^2/q.
$$
Тогда окончательно можем записать уравнение огибающей для круглого пучка радиуса $r$ с распределением Капчинского-Владимирского с внешней фокусировкой линейными полями:
$$
\displaystyle r'' + \frac{1}{\beta^2\gamma} \gamma' r' + \frac{1}{2\beta^2\gamma}\gamma''r + kr - \frac{2I}{I_a (\beta\gamma)^3}\frac{1}{r} - \frac{\epsilon^2}{r^3} = 0 ;
$$
и для эллиптического пучка:
$$
\begin{equation*}
 \begin{cases}
   \displaystyle a'' + \frac{1}{\beta^2\gamma} \gamma' a' + \frac{1}{2\beta^2\gamma}\gamma''a + k_xa - \frac{4I}{I_a (\beta\gamma)^3}\frac{1}{(a+b)} - \frac{\epsilon_x^2}{a^3} = 0 ,
   \\
   \displaystyle b'' + \frac{1}{\beta^2\gamma} \gamma' b' + \frac{1}{2\beta^2\gamma}\gamma''b + k_yb - \frac{4I}{I_a (\beta\gamma)^3}\frac{1}{(a+b)} - \frac{\epsilon_y^2}{b^3} = 0.
 \end{cases}
\end{equation*}
$$

### Магнитные линзы
Фокусирующими элементами могут являться: соленоиды и магнитные квадрупольные линзы.

### Соленоиды

$k_x = k_y = k_s$ - жесткость соленоида:
$$
k_s = \left ( \frac{eB_z}{2m_ec\beta\gamma} \right )^2 = \left ( \frac{e B_z}{2\beta\gamma\cdot 0.511\cdot 10^6 e \cdot \mathrm{volt}/c} \right )^2 =
\left ( \frac{cB_z[\mathrm{T}]}{2\beta\gamma\cdot 0.511\cdot 10^6 \cdot \mathrm{volt}} \right )^2.
$$

### Квадруполи

$k_q = \displaystyle\frac{eG}{pc}$ - жесткость квадруполя, где $G = \displaystyle\frac{\partial B_x}{\partial y} = \displaystyle\frac{\partial B_y}{\partial x}$ - градиент магнитного поля, причем $k_x = k_q, k_y = -k_q.$
$$
k_q = \left ( \frac{eG}{m_ec\beta\gamma} \right ) = \left ( \frac{eG}{\beta\gamma\cdot 0.511\cdot 10^6 e \cdot \mathrm{volt}/c} \right ) =
\left ( \frac{cG}{\beta\gamma\cdot 0.511\cdot 10^6 \cdot \mathrm{volt}} \right ).
$$

### Продольная динамика пучка

Уравнение на продольную динамику пучка можно решить независимо от уравнения на огибающую, чтобы в уравнении на огибающую уже использовать готовую функцию энергии пучка от $z$. Считая, что скорость электрона достаточно близка к скорости света и следовательно его продольная координата $z \approx ct$, а импульс $p_z \approx \gamma mc$

$$
\frac{d\gamma}{dz} \approx \frac{eE_z}{mc^2},
$$

Тогда достаточно один раз проинтегрировать функцию $E_z(z)$:

## Решение уравнения огибающей для эллиптического пучка  с фокусирующими элементами

Уравнение огибающей для эллиптического пучка с полуосями $a, b$ с распределением Капчинского-Владимирского с внешней фокусировкой линейными полями:
$$
\begin{equation*}
 \begin{cases}
   \displaystyle a'' + \frac{1}{\beta^2\gamma} \gamma' a' + \frac{1}{2\beta^2\gamma}\gamma''a + k_qa - \frac{2P}{(a+b)} - \frac{\epsilon_x^2}{a^3} = 0 ,
   \\
   \displaystyle b'' + \frac{1}{\beta^2\gamma} \gamma' b' + \frac{1}{2\beta^2\gamma}\gamma''b - k_qb - \frac{2P}{(a+b)} - \frac{\epsilon_y^2}{b^3} = 0.
 \end{cases}
\end{equation*}
$$


Пусть $\displaystyle x = \frac{da}{dz},  y = \frac{db}{dz},  \displaystyle \frac{d\gamma }{dz}\approx \frac{e E_z}{m c^2}, k = \frac{1}{R} - $ кривизна, тогда
$$
\displaystyle
\left\{\begin{matrix}
\displaystyle \frac{dx}{dz}= -  \frac{1}{\beta^2\gamma} \gamma' a' - \frac{1}{2\beta^2\gamma}\gamma''a - k_qa + \frac{2P}{(a+b)} + \frac{\epsilon_x^2}{a^3}\\ 
 \displaystyle\frac{da}{dz} =  x\\
 \displaystyle \frac{dy}{dz}= -  \frac{1}{\beta^2\gamma} \gamma' b' - \frac{1}{2\beta^2\gamma}\gamma''b + k_qb + \frac{2P}{(a+b)} + \frac{\epsilon_x^2}{b^3}\\ 
 \displaystyle\frac{db}{dz} =  y
\end{matrix}\right.
$$
Пусть $\vec X =
\begin{bmatrix}
x \\
a \\
y \\
b
\end{bmatrix} $, теперь составим дифференциальное уравнение $X' = F(X).$



### Литература

1. Дж. Лоусон "Физика пучков заряженных частиц"
2. Н.А. Винокуров "Лекции по электронной оптике для ускорительных физиков" 
