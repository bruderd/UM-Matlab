function vf2 = vf_sysid(in1,u1)
%VF_SYSID
%    VF2 = VF_SYSID(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    03-Aug-2018 13:23:39

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = u1.^2;
t3 = x1.^2;
t4 = x2.^2;
t5 = x3.^2;
t6 = x4.^2;
vf2 = [t2.*(-5.443405424768447e-6)+t3.*1.141049730879414e-4+t4.*1.141865464948433e-4-t5.*2.36562599899533e-4+t6.*1.654714136572074e-4-u1.*1.495958786174523e-4+x1.*1.668059768529475e-3-x2.*1.668343797543184e-3+x3.*1.002891697728362-x4.*2.621513814879435e-3-u1.*x1.*1.223011746188273e-4+u1.*x2.*1.223143949510788e-4+u1.*x3.*1.615407886146018e-4-u1.*x4.*1.496504581433841e-4-x1.*x2.*2.282913971527393e-4-x1.*x3.*8.651598598810277e-6-x1.*x4.*1.209187891882972e-4+x2.*x3.*8.527666341037952e-6+x2.*x4.*1.210151315486005e-4+x3.*x4.*6.088504563012113e-5+1.388377810452884e-4;t2.*2.215312213519521e-5-t3.*7.346684842165702e-5-t4.*7.358574023088793e-5+t5.*5.054926645727849e-4-t6.*9.259130835743919e-5+u1.*6.069478759492076e-5-x1.*2.082725591306488e-3+x2.*2.082732599410169e-3-x3.*3.550480945211133e-3+x4.*1.003488889657681+u1.*x1.*6.377883684540486e-5-u1.*x2.*6.380115850487705e-5-u1.*x3.*3.08568266070908e-4+u1.*x4.*1.348560663812519e-4+x1.*x2.*1.470525198418605e-4+x1.*x3.*7.02779549524668e-4-x1.*x4.*2.265494015305093e-4-x2.*x3.*7.028046480801396e-4+x2.*x4.*2.266096004613093e-4-x3.*x4.*2.47081636570145e-4-1.111866168230911e-3;t2.*(-3.159149582165911e-3)-t3.*1.095710800383176e-1-t4.*1.098475903050407e-1+t5.*6.280254893947361e-2-t6.*6.671846093885103e-2+u1.*2.385490877807085e-1-x1.*1.018666351310848e1+x2.*1.018675469597968e1-x3.*3.48868057243742e-1+x4.*1.203071786008849e-1-u1.*x1.*3.813144368172088e-1+u1.*x2.*3.813178626565301e-1-u1.*x3.*6.871125903577979e-2+u1.*x4.*6.48273402775898e-2+x1.*x2.*2.194185236389452e-1+x1.*x3.*1.56520275594599-x1.*x4.*2.891088756140378e-2-x2.*x3.*1.565239686304147+x2.*x4.*2.893502475151508e-2+x3.*x4.*6.439658157591067e-4+7.321036828360075e-1;t2.*3.779185717073789e-3+t3.*6.200802860925454e-1+t4.*6.218131341392521e-1-t5.*1.644038596328855e-1+t6.*7.929641026598239e-2-u1.*2.028166395922796e-1+x1.*1.902569821949058e1-x2.*1.902599319724938e1+x3.*8.865859071730644e-1-x4.*7.480618331184643e-1+u1.*x1.*4.953400918144832e-1-u1.*x2.*4.953311413957109e-1+u1.*x3.*1.031944003546304e-1-u1.*x4.*6.69155257514991e-2-x1.*x2.*1.241893205083639-x1.*x3.*7.073488825750539e-1-x1.*x4.*1.656023870220357+x2.*x3.*7.073734684382041e-1+x2.*x4.*1.656001164342728+x3.*x4.*6.476399172974634e-2-9.895163615655231e-1];
