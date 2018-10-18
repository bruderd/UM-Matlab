function vf2 = vf_koopman(in1,in2)
%VF_KOOPMAN
%    VF2 = VF_KOOPMAN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    18-Oct-2018 12:47:12

u1 = in2(1,:);
u2 = in2(2,:);
u3 = in2(3,:);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
t2 = u1.^2;
t3 = u2.^2;
t4 = u3.^2;
t5 = x1.^2;
t6 = x2.^2;
t7 = x3.^2;
t8 = t2.^2;
t9 = t3.^2;
t10 = t4.^2;
t11 = t5.^2;
t12 = t6.^2;
t13 = t7.^2;
vf2 = [t2.*(-8.359964493789223e-4)-t3.*1.441942748613769e-1+t4.*4.372690968781254e-2-t5.*7.149305373664133e-2+t6.*1.748208527829034e-2-t7.*2.817818177356612e-2-t8.*1.035128308739788e-1+t9.*2.136318118441143e-1-t10.*7.224390515600053e-2-t11.*6.864380389781903e-2-t12.*2.895732050602355e-3-t13.*2.267246777523467e-2-u1.*5.305952224744722e-3+u2.*9.727141856020356e-2-u3.*3.722102093495366e-2-x1.*1.142087478336027e-1+x2.*2.159822298124793e-2-x3.*6.2721453005044e-2+t2.*t3.*1.011112288184067e-1+t2.*t4.*1.068911821093819e-1-t2.*t5.*1.177222771797188e-1-t3.*t4.*1.622049981794422e-2-t2.*t6.*3.572009687377608e-2+t3.*t5.*1.093911189002047e-1+t2.*t7.*6.59777112171764e-2+t3.*t6.*1.111242353147264e-1+t4.*t5.*4.60176276263878e-2+t3.*t7.*1.441237793314253e-2+t4.*t6.*5.426080778566454e-2+t4.*t7.*3.363392576250435e-2+t5.*t6.*4.274229075961524e-2+t5.*t7.*7.414316640445058e-2+t6.*t7.*3.528556538533605e-2-t2.*u1.*2.556824671352594e-2+t2.*u2.*8.673232677413335e-2+t3.*u1.*9.76661843579992e-3+t2.*u3.*9.859851660936197e-3-t3.*u2.*1.094033967870036e-1-t4.*u1.*6.915962226250207e-2-t3.*u3.*3.365742885237513e-3+t4.*u2.*2.960518120053011e-2+t5.*u1.*1.824808705697311e-1+t4.*u3.*2.874794965027947e-2-t5.*u2.*9.278355289287896e-2+t6.*u1.*3.363893628384188e-4+t5.*u3.*1.457164723837869e-1-t6.*u2.*1.693911855475646e-1-t7.*u1.*3.932296704723856e-3+t6.*u3.*2.223120447103876e-2+t7.*u2.*4.754550075378456e-2-t7.*u3.*5.819594520281832e-2+u1.*u2.*2.350650853006052e-2-u1.*u3.*8.532078392316514e-2+u2.*u3.*1.436466278517778e-1-t2.*x1.*3.039222906866657e-2-t2.*x2.*1.73596587323655e-3-t3.*x1.*1.80428231321914e-2+t2.*x3.*7.507819263143341e-2-t3.*x2.*1.395094184831605e-2+t4.*x1.*2.22901025674935e-2+t3.*x3.*6.149907958081541e-2-t4.*x2.*2.602274481193903e-2-t5.*x1.*2.467497760678414e-1+t4.*x3.*2.239040171522261e-2+t5.*x2.*6.025226328930242e-2-t6.*x1.*3.020155080745301e-1-t5.*x3.*1.630992411724531e-1-t6.*x2.*5.952540339478601e-2+t7.*x1.*5.22113699541733e-2+t6.*x3.*6.526800501872398e-2+t7.*x2.*9.843967681519951e-2+t7.*x3.*2.538634321780163e-2+u1.*x1.*5.633226621694185e-2+u1.*x2.*6.809554504006365e-2+u2.*x1.*2.032215525133033e-1+u1.*x3.*1.945722411250377e-1-u2.*x2.*9.091178392310754e-2-u3.*x1.*8.049436781785001e-2-u2.*x3.*2.65909521599903e-1-u3.*x2.*1.055129826256679e-2+u3.*x3.*1.42710558231934e-1+x1.*x2.*5.122031795945214e-2-x1.*x3.*1.951576660727992e-1+x2.*x3.*1.882397514695546e-2+t2.*u1.*u2.*3.317207748140112e-4+t2.*u1.*u3.*4.183864274632554e-2+t3.*u1.*u2.*1.596463752901696e-3+t2.*u2.*u3.*5.394590021834992e-3-t3.*u1.*u3.*1.561598112563563e-1+t4.*u1.*u2.*3.025247083413285e-3+t3.*u2.*u3.*9.3482256379975e-2+t4.*u1.*u3.*1.666416743671203e-2+t5.*u1.*u2.*6.286299872052883e-3-t4.*u2.*u3.*1.238121156872529e-1-t5.*u1.*u3.*2.785821886576586e-1+t6.*u1.*u2.*6.89614266814256e-2-t5.*u2.*u3.*1.53982319183206e-1+t6.*u1.*u3.*5.327202454115709e-2-t7.*u1.*u2.*4.318600713097151e-2-t6.*u2.*u3.*1.236670071955861e-1+t7.*u1.*u3.*7.884276008307256e-2+t7.*u2.*u3.*6.547779476734053e-2-u1.*u2.*u3.*3.60903425119288e-2-t2.*u1.*x1.*1.048705720986048e-1-t2.*u1.*x2.*8.01841348245567e-2-t2.*u2.*x1.*4.48859913128201e-2-t3.*u1.*x1.*6.095021245733503e-2-t2.*u1.*x3.*2.031731975578576e-2-t2.*u2.*x2.*2.072347004962755e-1+t2.*u3.*x1.*1.392614767696791e-1+t3.*u1.*x2.*1.665738042541133e-1+t3.*u2.*x1.*1.661478485393195e-1+t4.*u1.*x1.*1.179310880695809e-2+t2.*u2.*x3.*5.971739836556119e-3+t2.*u3.*x2.*3.756741420953291e-2+t3.*u1.*x3.*1.174082615376257e-1+t3.*u2.*x2.*9.044965723344652e-3-t3.*u3.*x1.*1.636326596550252e-1+t4.*u1.*x2.*8.748422619557104e-2-t4.*u2.*x1.*1.409421178083573e-2+t5.*u1.*x1.*2.952721706156149e-1-t2.*u3.*x3.*5.112415775077703e-2+t3.*u2.*x3.*8.45351135084399e-2+t3.*u3.*x2.*1.72966974494606e-1-t4.*u1.*x3.*1.045820060223618e-1+t4.*u2.*x2.*2.502959838949417e-2-t4.*u3.*x1.*6.60923965220466e-2+t5.*u1.*x2.*1.460695764287183e-1+t5.*u2.*x1.*3.783484073414949e-1+t6.*u1.*x1.*2.103977226124094e-1-t3.*u3.*x3.*9.335551377351925e-2-t4.*u2.*x3.*1.328997305904966e-1-t4.*u3.*x2.*9.650906784695731e-3-t5.*u1.*x3.*7.410680748646825e-2+t5.*u2.*x2.*9.698574617058585e-3-t5.*u3.*x1.*2.055921890180705e-1+t6.*u1.*x2.*5.854187918251548e-2+t6.*u2.*x1.*3.836680377513569e-1-t7.*u1.*x1.*5.289269179374079e-2-t4.*u3.*x3.*4.237037065498558e-2+t5.*u2.*x3.*2.180772827552745e-1-t5.*u3.*x2.*2.034737428626037e-1-t6.*u1.*x3.*7.893411598347828e-2+t6.*u2.*x2.*3.543993049250764e-2+t6.*u3.*x1.*1.287885059133864e-1-t7.*u1.*x2.*2.668609577674152e-2-t7.*u2.*x1.*6.306585969600982e-2+t5.*u3.*x3.*2.276623147131451e-2+t6.*u2.*x3.*1.169489499691694e-1+t6.*u3.*x2.*9.621795765054132e-2-t7.*u1.*x3.*1.16988912739291e-1-t7.*u2.*x2.*2.013332712798933e-2+t7.*u3.*x1.*8.420207021485691e-2-t6.*u3.*x3.*1.161537355405916e-1+t7.*u2.*x3.*5.941887026886366e-2-t7.*u3.*x2.*3.112039576123879e-2-t7.*u3.*x3.*5.208035624188644e-2+u1.*u2.*x1.*1.731647262678471e-1+u1.*u2.*x2.*8.303681181471132e-3+u1.*u3.*x1.*3.878043159758484e-2+u1.*u2.*x3.*7.446326631970315e-2-u1.*u3.*x2.*5.709110463811544e-2+u2.*u3.*x1.*7.303471204317956e-2+u1.*u3.*x3.*6.782506146151118e-2-u2.*u3.*x2.*2.379630058605039e-2-u2.*u3.*x3.*8.354612814135347e-2-t2.*x1.*x2.*4.407981057677821e-2+t2.*x1.*x3.*5.504999077051885e-2-t3.*x1.*x2.*1.895494352462412e-1+t2.*x2.*x3.*3.116582386176378e-2+t3.*x1.*x3.*4.800955016917561e-1-t4.*x1.*x2.*8.752396581445296e-2+t3.*x2.*x3.*8.779408325310994e-2-t4.*x1.*x3.*2.069356417429233e-1+t5.*x1.*x2.*2.258564572297257e-1-t4.*x2.*x3.*9.614725208759172e-3-t5.*x1.*x3.*2.867320696169976e-2+t6.*x1.*x2.*1.032343240477334e-1+t5.*x2.*x3.*5.020479037237225e-2+t6.*x1.*x3.*2.210289948593492e-2-t7.*x1.*x2.*3.172088903765187e-2+t6.*x2.*x3.*4.022207263793562e-2-t7.*x1.*x3.*1.098484249668708e-1-t7.*x2.*x3.*2.12994803655485e-2+u1.*x1.*x2.*3.537530025352604e-2+u1.*x1.*x3.*1.463691411870554e-1+u2.*x1.*x2.*1.320631098950495e-1-u1.*x2.*x3.*9.816163005951348e-3+u2.*x1.*x3.*2.789726046648249e-1-u3.*x1.*x2.*1.038710259635793e-1+u2.*x2.*x3.*4.382609555463295e-2+u3.*x1.*x3.*4.519418338729657e-2+u3.*x2.*x3.*8.630259974591779e-2+x1.*x2.*x3.*1.495211814498911e-1-u1.*u2.*u3.*x1.*2.287036439383466e-1+u1.*u2.*u3.*x2.*2.362605576807458e-1+u1.*u2.*u3.*x3.*1.233356351102768e-1+u1.*u2.*x1.*x2.*3.588262646619979e-2+u1.*u2.*x1.*x3.*6.923709061579213e-2-u1.*u3.*x1.*x2.*6.998194329081471e-3-u1.*u2.*x2.*x3.*5.0487442822541e-2+u1.*u3.*x1.*x3.*9.88821003791335e-2+u2.*u3.*x1.*x2.*2.018000624199925e-1-u1.*u3.*x2.*x3.*6.010047175899635e-3-u2.*u3.*x1.*x3.*6.078340729688164e-2+u2.*u3.*x2.*x3.*8.358469085399277e-2-u1.*x1.*x2.*x3.*3.335096245109052e-2+u2.*x1.*x2.*x3.*5.143725729955766e-2+u3.*x1.*x2.*x3.*8.759795159488498e-2-1.290065189039369e-2;t2.*(-3.577090253146072e-2)+t3.*6.303704703727467e-2-t4.*1.298385466489502e-1+t5.*2.184151257232681e-3-t6.*9.929290364113182e-2+t7.*7.901299672182245e-2+t8.*4.124117076450102e-3-t9.*2.457556872417646e-1+t10.*3.434335165016603e-1+t11.*1.192983750616609e-1+t12.*1.935171313746999e-2-t13.*1.247614256899619e-2-u1.*1.988105499944411e-2-u2.*9.116083339785383e-2+u3.*6.567308300184865e-2-x1.*1.395111233652807e-2-x2.*1.103344885847209e-1+x3.*2.761521413275184e-2+t2.*t3.*8.602435867281333e-2-t2.*t4.*6.693870148245043e-2+t2.*t5.*8.857991500467084e-3+t3.*t4.*4.87406358993922e-2-t2.*t6.*1.560740429597287e-1-t3.*t5.*3.16199082793901e-2-t2.*t7.*2.853674319277521e-2-t3.*t6.*8.198196328066634e-2+t4.*t5.*1.489528100501541e-1+t3.*t7.*1.112564877603138e-1+t4.*t6.*1.569650627535155e-1-t4.*t7.*1.01637587788695e-1+t5.*t6.*1.278145729968196e-1-t5.*t7.*7.459197723407401e-2+t6.*t7.*3.530676212220673e-2+t2.*u1.*1.923305682306585e-2-t2.*u2.*5.62199960401738e-2+t3.*u1.*1.821531823064357e-1-t2.*u3.*3.264219515573336e-2+t3.*u2.*5.710037389214305e-3-t4.*u1.*4.919992324058965e-2+t3.*u3.*6.088176775838908e-2+t4.*u2.*1.08549107231387e-1+t5.*u1.*2.215467646916504e-2-t4.*u3.*1.130664085053472e-1-t5.*u2.*1.337466543295343e-2+t6.*u1.*5.842198652421696e-2+t5.*u3.*7.205622353950129e-4+t6.*u2.*3.037082620112429e-2-t7.*u1.*4.441556537189268e-2-t6.*u3.*1.383748162164368e-1+t7.*u2.*9.200030511928241e-2+t7.*u3.*6.098243098333597e-2-u1.*u2.*1.119399432163338e-1+u1.*u3.*1.367537151768251e-1+u2.*u3.*2.025421199555428e-2+t2.*x1.*4.615716379749594e-2-t2.*x2.*1.30841356706879e-1+t3.*x1.*5.57985895843292e-2+t2.*x3.*1.106433697454164e-2-t3.*x2.*1.123673707131528e-1+t4.*x1.*9.721600821647332e-2-t3.*x3.*7.787005744958687e-2-t4.*x2.*1.544028316562157e-1+t5.*x1.*5.144323624259193e-3+t4.*x3.*5.876015660163121e-2-t5.*x2.*1.221780995441866e-1+t6.*x1.*9.853591667487134e-2+t5.*x3.*1.733852275232567e-1-t6.*x2.*1.283562364051976e-1-t7.*x1.*8.285494312207782e-2-t6.*x3.*5.274056093595827e-2+t7.*x2.*2.238963486980943e-2+t7.*x3.*1.418037050101916e-3+u1.*x1.*1.021854748035612e-1-u1.*x2.*5.209476565251826e-2-u2.*x1.*1.093770447399747e-1+u1.*x3.*1.683342496214878e-1-u2.*x2.*2.765430666267029e-2+u3.*x1.*6.064957907461829e-2+u2.*x3.*3.466074672715128e-2+u3.*x2.*3.410821199505944e-1-u3.*x3.*2.671101019593785e-1+x1.*x2.*5.945507764344199e-3+x1.*x3.*5.520168489148478e-2-x2.*x3.*8.049358109531546e-2-t2.*u1.*u2.*8.645884737315052e-2-t2.*u1.*u3.*2.798244950117738e-2+t3.*u1.*u2.*2.678206831173278e-2-t2.*u2.*u3.*3.055364287516271e-2-t3.*u1.*u3.*6.73020971100992e-3-t4.*u1.*u2.*9.744061020386185e-2+t3.*u2.*u3.*1.524576406893334e-1+t4.*u1.*u3.*3.905658024441665e-2+t5.*u1.*u2.*1.141523970414722e-1+t4.*u2.*u3.*9.374163725935154e-3-t5.*u1.*u3.*2.818054104874732e-2-t6.*u1.*u2.*2.527389953112041e-2+t5.*u2.*u3.*8.253899363270381e-3+t6.*u1.*u3.*7.393340046216858e-2-t7.*u1.*u2.*3.42971842075842e-2+t6.*u2.*u3.*1.151075338873534e-1-t7.*u1.*u3.*1.74165599020031e-2-t7.*u2.*u3.*1.117569217768195e-1-u1.*u2.*u3.*8.009667915767406e-2+t2.*u1.*x1.*3.055257033484902e-2+t2.*u1.*x2.*3.444870133868498e-1-t2.*u2.*x1.*1.297090201929545e-1-t3.*u1.*x1.*6.591419525712872e-2-t2.*u1.*x3.*1.171191423966237e-1+t2.*u2.*x2.*9.59066263866856e-2+t2.*u3.*x1.*3.624380891006301e-2+t3.*u1.*x2.*8.01768565292124e-2+t3.*u2.*x1.*1.891936638129181e-1-t4.*u1.*x1.*1.679305078679233e-2-t2.*u2.*x3.*8.563073126851284e-2-t2.*u3.*x2.*1.119485596626997e-1+t3.*u1.*x3.*8.699421912995173e-2+t3.*u2.*x2.*4.082554987253703e-2+t3.*u3.*x1.*1.264701241370571e-1-t4.*u1.*x2.*8.995501798531816e-2-t4.*u2.*x1.*1.373005014680936e-1+t5.*u1.*x1.*8.269118284171248e-2-t2.*u3.*x3.*9.597371720776882e-2+t3.*u2.*x3.*3.564384590539828e-2-t3.*u3.*x2.*9.630419960764792e-3-t4.*u1.*x3.*7.23230348086083e-3-t4.*u2.*x2.*7.567534457080044e-2-t4.*u3.*x1.*1.162750003988837e-1+t5.*u1.*x2.*8.743696007301525e-2-t5.*u2.*x1.*2.032810764608991e-1+t6.*u1.*x1.*1.489683340407998e-1+t3.*u3.*x3.*1.325630670280918e-1+t4.*u2.*x3.*2.503144400103164e-1+t4.*u3.*x2.*2.64256372765746e-2-t5.*u1.*x3.*4.516796731936366e-2+t5.*u2.*x2.*1.882699829404903e-1-t5.*u3.*x1.*5.642440143377474e-2+t6.*u1.*x2.*4.89544298181575e-2-t6.*u2.*x1.*2.765509241850938e-1-t7.*u1.*x1.*7.11303450729585e-2+t4.*u3.*x3.*1.190080603109556e-1+t5.*u2.*x3.*1.731054582895799e-2+t5.*u3.*x2.*3.522371055622573e-1-t6.*u1.*x3.*2.498849862648761e-1+t6.*u2.*x2.*1.069792696841774e-1-t6.*u3.*x1.*5.892541226928083e-2-t7.*u1.*x2.*3.548401350821761e-2+t7.*u2.*x1.*2.303388882946247e-1+t5.*u3.*x3.*4.415664274766496e-2-t6.*u2.*x3.*4.323447460015676e-2+t6.*u3.*x2.*2.196524229462196e-1-t7.*u1.*x3.*5.529208626923001e-3+t7.*u2.*x2.*2.224269869311038e-3+t7.*u3.*x1.*8.717374322891393e-2+t6.*u3.*x3.*2.01622513248766e-1+t7.*u2.*x3.*3.82390742818114e-2-t7.*u3.*x2.*5.420006427838473e-2+t7.*u3.*x3.*2.566509364952758e-2-u1.*u2.*x1.*2.652517653452116e-2-u1.*u2.*x2.*7.364741786773304e-2-u1.*u3.*x1.*1.087513561705529e-2+u1.*u2.*x3.*2.741195293544978e-2+u1.*u3.*x2.*1.666418317142428e-1-u2.*u3.*x1.*1.532304881557061e-1-u1.*u3.*x3.*1.791666175307282e-1+u2.*u3.*x2.*1.435239418254583e-1+u2.*u3.*x3.*7.587676464161507e-2-t2.*x1.*x2.*1.522058741715267e-2+t2.*x1.*x3.*5.184300731811525e-2+t3.*x1.*x2.*1.743310420716534e-1+t2.*x2.*x3.*2.692790915277407e-1-t3.*x1.*x3.*6.775201309488835e-2-t4.*x1.*x2.*5.052481147168034e-2-t3.*x2.*x3.*6.077283776011207e-3-t4.*x1.*x3.*7.091926173234647e-2-t5.*x1.*x2.*1.210376612159424e-1+t4.*x2.*x3.*3.485836226119366e-1-t5.*x1.*x3.*2.610000923335246e-2-t6.*x1.*x2.*5.519573852930998e-3+t5.*x2.*x3.*1.264661663552758e-1+t6.*x1.*x3.*8.659338830541496e-2+t7.*x1.*x2.*3.920815375189608e-2-t6.*x2.*x3.*1.04912321610586e-1-t7.*x1.*x3.*2.004992200467157e-2-t7.*x2.*x3.*1.865796520819145e-1+u1.*x1.*x2.*1.133232400329377e-2+u1.*x1.*x3.*6.369209257783999e-2-u2.*x1.*x2.*1.267170200731192e-1+u1.*x2.*x3.*9.19779869787071e-2-u2.*x1.*x3.*1.457863934974786e-1+u3.*x1.*x2.*5.135844030227468e-2+u2.*x2.*x3.*6.975821419659138e-2-u3.*x1.*x3.*8.862626567905629e-2+u3.*x2.*x3.*2.153068033336624e-1-x1.*x2.*x3.*1.543305213764005e-2-u1.*u2.*u3.*x1.*9.824402118301451e-2-u1.*u2.*u3.*x2.*1.113586238806629e-1+u1.*u2.*u3.*x3.*7.900880442184496e-3+u1.*u2.*x1.*x2.*9.379710917982143e-2+u1.*u2.*x1.*x3.*9.682006830380441e-2-u1.*u3.*x1.*x2.*4.748310781344807e-3-u1.*u2.*x2.*x3.*3.934035307873856e-2+u1.*u3.*x1.*x3.*1.626918804083099e-2-u2.*u3.*x1.*x2.*7.177251061158658e-2+u1.*u3.*x2.*x3.*8.393013642762705e-2+u2.*u3.*x1.*x3.*4.209374148168808e-2+u2.*u3.*x2.*x3.*5.362700839610121e-2+u1.*x1.*x2.*x3.*5.796461515643873e-3+u2.*x1.*x2.*x3.*8.038784313441291e-2-u3.*x1.*x2.*x3.*1.617675005179036e-2+1.666390901306916e-2;t2.*3.97874445828201e-2+t3.*3.709417550132797e-2-t4.*2.778054829591857e-3+t5.*1.8814544265521e-1+t6.*1.202637717784138e-1-t7.*3.049218156075983e-1-t8.*5.542200228425129e-2-t9.*5.452961053031193e-2-t10.*8.702225280062141e-2+t11.*4.767026814161855e-2-t12.*1.354226591967382e-1-t13.*1.65794703430181e-1+u1.*2.504922190555389e-2-u2.*3.406062892832796e-2+u3.*2.228103478448874e-2+x1.*1.774024124933661e-2-x2.*1.699536578660539e-2-x3.*5.384116789322031e-2-t2.*t3.*9.74048841320963e-2-t2.*t4.*1.08002471570225e-1+t2.*t5.*1.008931626445446e-1-t3.*t4.*6.056433890916848e-2-t2.*t6.*5.037030413485886e-2+t3.*t5.*9.305912413532602e-2+t2.*t7.*1.144628743612167e-1-t3.*t6.*1.47766689136418e-1-t4.*t5.*1.79653589828595e-1+t3.*t7.*2.328535274754562e-1+t4.*t6.*2.130216730561524e-1+t4.*t7.*9.217615563940355e-2-t5.*t6.*9.118446374125853e-2+t5.*t7.*2.904307455803489e-2-t6.*t7.*6.777894342885063e-2-t2.*u1.*9.715755367102671e-3+t2.*u2.*1.305844421434877e-2+t3.*u1.*1.548222462404729e-2-t2.*u3.*6.209182475119213e-2+t3.*u2.*2.087955121325785e-2-t4.*u1.*4.058355054480347e-2-t3.*u3.*1.788643177964264e-3+t4.*u2.*3.637719130083585e-2-t5.*u1.*7.732371178553303e-2+t4.*u3.*7.659503296075437e-2-t5.*u2.*2.635676251067178e-1+t6.*u1.*9.204288400233841e-2+t5.*u3.*1.392682419457996e-1+t6.*u2.*1.406638556458716e-1+t7.*u1.*1.406341970806952e-2-t6.*u3.*5.234199324193292e-2+t7.*u2.*2.472424467145631e-1+t7.*u3.*1.340273924966263e-1+u1.*u2.*6.405662257692962e-3-u1.*u3.*6.85846938985704e-2-u2.*u3.*6.000575602194617e-2+t2.*x1.*3.263076213182007e-2+t2.*x2.*1.091882741771197e-1+t3.*x1.*1.400601219799707e-2-t2.*x3.*4.268223534085094e-2+t3.*x2.*4.12614168441514e-2-t4.*x1.*7.084113710457461e-2-t3.*x3.*1.006629530323127e-1-t4.*x2.*2.146732857339806e-2-t5.*x1.*4.985098022452072e-2-t4.*x3.*3.592249942445097e-2-t5.*x2.*1.124881477479895e-2-t6.*x1.*1.123114134033753e-2+t5.*x3.*2.078330150703873e-2+t6.*x2.*3.620344208767937e-2-t7.*x1.*4.910307158670896e-2-t6.*x3.*1.631454446979038e-1-t7.*x2.*7.846002793023122e-3-t7.*x3.*2.365453109749096e-1+u1.*x1.*4.115357654980056e-2+u1.*x2.*9.059861665184483e-2-u2.*x1.*6.736254235395468e-2+u1.*x3.*1.277425332272537e-1+u2.*x2.*8.333562779464843e-2+u3.*x1.*6.46440949233051e-3-u2.*x3.*5.746785979315811e-3-u3.*x2.*8.901479150693524e-2+u3.*x3.*8.793235098836655e-2+x1.*x2.*3.581619731560602e-3+x1.*x3.*5.610513271024976e-2+x2.*x3.*5.371592791009537e-2+t2.*u1.*u2.*2.826071951722038e-2+t2.*u1.*u3.*5.909244855145328e-2-t3.*u1.*u2.*4.901455141275212e-2-t2.*u2.*u3.*3.843760374253478e-2+t3.*u1.*u3.*1.594591265906898e-1-t4.*u1.*u2.*1.432780690330097e-2-t3.*u2.*u3.*5.555797598461476e-3+t4.*u1.*u3.*8.56193870023467e-3-t5.*u1.*u2.*1.461213573107195e-1+t4.*u2.*u3.*7.240192478577061e-2-t5.*u1.*u3.*2.723939561976899e-2-t6.*u1.*u2.*2.389980032607555e-2+t5.*u2.*u3.*1.124521062005691e-1-t6.*u1.*u3.*3.122901022738759e-1-t7.*u1.*u2.*1.538337330289408e-1-t6.*u2.*u3.*3.224842305762152e-1+t7.*u1.*u3.*2.562639485131625e-1-t7.*u2.*u3.*5.280503887032365e-2+u1.*u2.*u3.*5.973518802549466e-3+t2.*u1.*x1.*1.027813767678828e-1-t2.*u1.*x2.*1.52178558969768e-2-t2.*u2.*x1.*1.732954628783171e-1+t3.*u1.*x1.*4.984471702426411e-3+t2.*u1.*x3.*9.936278523391366e-2+t2.*u2.*x2.*6.84840725428858e-2+t2.*u3.*x1.*2.826462857226039e-2-t3.*u1.*x2.*2.465710214715356e-2-t3.*u2.*x1.*3.110805704217557e-2-t4.*u1.*x1.*2.048786292562692e-2+t2.*u2.*x3.*1.622702154900291e-2+t2.*u3.*x2.*1.749381362443441e-2+t3.*u1.*x3.*1.73598399064486e-2+t3.*u2.*x2.*1.938898926758348e-1-t3.*u3.*x1.*1.099403348515113e-1+t4.*u1.*x2.*1.307615622761435e-1+t4.*u2.*x1.*2.172511638995953e-2-t5.*u1.*x1.*5.159370199551188e-3-t2.*u3.*x3.*1.652753465899112e-1+t3.*u2.*x3.*1.082915322568384e-1-t3.*u3.*x2.*2.749981304482253e-1-t4.*u1.*x3.*9.706536220428941e-2-t4.*u2.*x2.*2.77378291217942e-2+t4.*u3.*x1.*3.158949035931667e-2-t5.*u1.*x2.*1.579447374358361e-1+t5.*u2.*x1.*4.362846771969298e-2+t6.*u1.*x1.*1.196548299805241e-1+t3.*u3.*x3.*2.288671864852852e-3+t4.*u2.*x3.*1.308807624790269e-3-t4.*u3.*x2.*1.638654381330533e-1+t5.*u1.*x3.*2.458845281876162e-2-t5.*u2.*x2.*1.043343461149232e-1+t5.*u3.*x1.*2.210536770901309e-1-t6.*u1.*x2.*4.351260534037461e-2+t6.*u2.*x1.*9.598326646225536e-2-t7.*u1.*x1.*7.312746187746248e-2+t4.*u3.*x3.*2.710192740942046e-2+t5.*u2.*x3.*4.608139641186139e-2+t5.*u3.*x2.*1.292422630097553e-1-t6.*u1.*x3.*4.462064221102675e-2-t6.*u2.*x2.*1.374594742040732e-2-t6.*u3.*x1.*5.03561836088528e-3-t7.*u1.*x2.*1.477283591744344e-1+t7.*u2.*x1.*1.159843347895302e-1+t5.*u3.*x3.*1.340010298028544e-1+t6.*u2.*x3.*2.115350067892093e-2+t6.*u3.*x2.*3.544943894914231e-2+t7.*u1.*x3.*4.13880319210449e-3+t7.*u2.*x2.*1.1990004854686e-1-t7.*u3.*x1.*6.337036218641488e-2-t6.*u3.*x3.*1.416971171515206e-2+t7.*u2.*x3.*7.761115615789367e-2-t7.*u3.*x2.*5.679777464425247e-2+t7.*u3.*x3.*3.85091521327171e-2-u1.*u2.*x1.*5.493652070107986e-2-u1.*u2.*x2.*1.09647892931433e-1+u1.*u3.*x1.*8.116696390394378e-2-u1.*u2.*x3.*8.151148053689861e-2+u1.*u3.*x2.*7.812777589812773e-3-u2.*u3.*x1.*3.087413937022281e-2-u1.*u3.*x3.*6.27316780299748e-2-u2.*u3.*x2.*3.822458698770267e-2+u2.*u3.*x3.*1.894882573013943e-2+t2.*x1.*x2.*1.493708009419942e-1-t2.*x1.*x3.*7.853556690717528e-2-t3.*x1.*x2.*2.166888742541062e-1+t2.*x2.*x3.*4.511436883212502e-2+t3.*x1.*x3.*5.539222250797255e-2+t4.*x1.*x2.*1.254591201277778e-1-t3.*x2.*x3.*4.877822497229942e-2+t4.*x1.*x3.*2.417022174512083e-2-t5.*x1.*x2.*4.137225130306808e-3-t4.*x2.*x3.*4.244628093477463e-2-t5.*x1.*x3.*5.880827182903819e-3+t6.*x1.*x2.*1.028843502153554e-1-t5.*x2.*x3.*4.738708662620154e-2+t6.*x1.*x3.*1.611969177657812e-1-t7.*x1.*x2.*4.306548112435739e-2+t6.*x2.*x3.*1.006291761836872e-2+t7.*x1.*x3.*2.031068612163547e-2-t7.*x2.*x3.*3.715634918680299e-3-u1.*x1.*x2.*1.404539500065358e-1-u1.*x1.*x3.*6.657472455608814e-2+u2.*x1.*x2.*1.045852644072002e-1-u1.*x2.*x3.*3.979451313234737e-2+u2.*x1.*x3.*3.737981056075948e-2+u3.*x1.*x2.*7.280707448264676e-2+u2.*x2.*x3.*3.055646156897534e-2+u3.*x1.*x3.*4.204619639466858e-2+u3.*x2.*x3.*1.099532936629638e-2+x1.*x2.*x3.*1.379628184760572e-1+u1.*u2.*u3.*x1.*1.470604384239619e-1+u1.*u2.*u3.*x2.*3.754194200646711e-2-u1.*u2.*u3.*x3.*1.150999672578579e-1+u1.*u2.*x1.*x2.*1.060956502038704e-1-u1.*u2.*x1.*x3.*5.349493940840407e-2-u1.*u3.*x1.*x2.*2.529903037122223e-1-u1.*u2.*x2.*x3.*1.616958868637784e-2-u1.*u3.*x1.*x3.*1.425863142555242e-1+u2.*u3.*x1.*x2.*4.453896464172564e-2+u1.*u3.*x2.*x3.*1.00762479811381e-1+u2.*u3.*x1.*x3.*6.574824455636842e-2-u2.*u3.*x2.*x3.*1.784055226735543e-1+u1.*x1.*x2.*x3.*2.826361171229448e-3-u2.*x1.*x2.*x3.*6.64212459188314e-2-u3.*x1.*x2.*x3.*1.130888005599778e-1+2.145711356400931e-4];
