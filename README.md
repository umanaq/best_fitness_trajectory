# best_fitness_trajectory

运动学计算公式
【加加速 phase 1】
$$
v(t)=\int_0^t Jt\ \mathbb{dt}=\frac{Jt^2}{2}
\\ s(t)= \int_0^t \frac{Jt^2}{2}\ \mathbb{dt} = \frac{Jt^3}{6}
\\ 其中\ t\in[0, t_1]
$$

$$
\therefore v_1=v(t_1)=\frac{Jt_j^2}{2},\quad s_1=s(t_1)=\frac{Jt_j^3}{6}
$$

【匀加速 phase 2】
$$
v(t) = v_1+a_mt = v_1+Jt_1t
\\s(t) = \int_0^t [\frac{Jt_1^2}{2}+Jt_1t]\ \mathbb{dt}
=\frac{1}{2}Jt_1t(t+t_1)
\\ 其中\ t\in[0, t_2]
$$

$$
\therefore v_2 = v(t_2) = Jt_j(\frac{t_j}{2} + t_a),\quad s_2=s(t_2)=\frac{1}{2}Jt_jt_a(t_j+t_a)
$$

【减加速 phase 3】
$$
v(t)=v_2+\int_0^t (a_m-Jt)\ \mathbb{dt}=v_2+Jt_1t-\frac{Jt^2}{2}
\\s(t)=\int_0^t[(\frac{Jt_1^2}{2}+Jt_1t_2)+Jt_1t-\frac{Jt^2}{2}]\ \mathbb{dt}
\\=(\frac{Jt_1^2}{2}+Jt_1t_2)t+\frac{Jt_1t^2}{2}-\frac{Jt^3}{6}
\\= \frac{Jt_1t}{2} (t_1+2t_2+t) - \frac{Jt^3}{6}
\\ 其中\ t\in[0, t_3]
$$

$$
\therefore v_3=v(t_3) = Jt_j(t_j+t_a) ,\quad s_3=s(t_3)=Jt_j^2(\frac{5}{6}t_j+t_a)
$$

【匀速 phase 4】
$$
v(t)=v_3=Jt_1(t_1+t_2)
\\s(t)=\int_0^t Jt_1(t_1+t_2)\ \mathbb{dt}=Jt_1(t_1+t_2)t
\\ 其中\ t\in[0, t_4]
$$

$$
\therefore v_4=v(t_3) = Jt_j(t_j+t_a) ,\quad s_4=s(t_4)=Jt_jt_v(t_j+t_a)
$$

【加减速 phase 5】
$$
v(t) = v_4+\int_0^t (-\frac{J}{\alpha^2}t)\ \mathbb{dt}=v_4-\frac{Jt^2}{2\alpha^2}\\
s(t) = \int_0^t [Jt_1(t_1+t_2) -\frac{Jt^2}{2\alpha^2}] \ \mathbb{dt}\\
=Jt_1(t_1+t_2)t -\frac{Jt^3}{6\alpha^2}
\\ 其中\ t\in[0, t_5]
$$

$$
\therefore v_5 = v(t_5) = Jt_j(\frac{t_j}{2} + t_a) = v_2
\\s_5=s(t_5)=J\alpha t_j^2(t_j+t_a)-\frac{J\alpha t_j^3}{6}
=\alpha Jt_j^2(\frac{5}{6}t_j + t_a) = \alpha s_3
$$

【匀减速 phase 6】
$$
v(t) = v_5+\int_0^t(-\frac{Jt_5}{\alpha^2}) \ \mathbb{dt} =  v_5-\frac{Jt_5}{\alpha^2}t = v_5 - \frac{Jt_j}{\alpha}t\\
s(t) = \int_0^t [Jt_j(\frac{t_j}{2} + t_a)- \frac{Jt_j}{\alpha}t]\ \mathbb{dt}
= Jt_j(\frac{t_j}{2} + t_a)t - \frac{Jt_j}{2\alpha}t^2
\\=\frac{Jt_jt}{2}[t_j + 2t_a - \frac{t}{\alpha}]
\\ 其中\ t\in[0, t_6]
$$

$$
\therefore v_6=v(t_6)=v_2-Jt_jt_a=\frac{Jt_j^2}{2}=v_1
\\s_6=s(t_6)= Jt_j(\frac{t_j}{2} + t_a)\alpha t_a - \frac{Jt_j}{2}\alpha t_a^2
= \frac{\alpha}{2}Jt_jt_a(t_j+t_a)=\alpha s_2
$$

【减减速 phase 7】
$$
v(t)=v_6+\int_0^t (-\frac{Jt_5}{\alpha^2} +\frac{J}{\alpha^2}t) \ \mathbb{dt} = v_6-\frac{Jt_jt}{\alpha}+\frac{Jt^2}{2\alpha^2}
\\s(t)=\int_0^t (\frac{Jt_j^2}{2}-\frac{Jt_jt}{\alpha}+\frac{Jt^2}{2\alpha^2}) \ \mathbb{dt}
= \frac{Jt_j^2t}{2} - \frac{Jt_jt^2}{2\alpha} + \frac{Jt^3}{6\alpha^2}
\\=\frac{Jt}{2}(t_j^2-\frac{tt_j}{\alpha}+\frac{t^2}{3\alpha^2})
\\ 其中\ t\in[0, t_7]
$$

$$
\therefore v_7= v_1-Jt_j^2+\frac{Jt_j^2}{2}=0,\quad s_7=s(t_7)=\frac{J\alpha t_j^3}{2} - \frac{J\alpha^2t_j^3}{2\alpha} + \frac{J\alpha^3t_j^3}{6\alpha^2}=\frac{\alpha Jt_j^3}{6}=\alpha s_1
$$

