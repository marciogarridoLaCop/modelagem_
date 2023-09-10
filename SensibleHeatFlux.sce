clear
clc

Dados = csvRead('/Users/marciogarrido/resltados_experimento_artigo/DADOS_exportados40.csv')

// Parâmetros
Cp = 1013       // Calor específico do ar [J/kg oC]
Cpm = 29.30     // Calor específico do ar - molar [J/mol oC]
Ro = 1.146      // Densidade do ar [kg/m²]
Lc = 0.155      // Diâmetro da placa [m]
Vcf = 1.64E-5   // Coeficiente de viscosidade cinemática do ar [m²/s]
Dh = 2.67E-2    // Difusividade térmica [W/m K]

nrows = size (Dados,"r")
j = 0

// Filtrar os valores negativos do vento
for i = 2:nrows
    if Dados(i,8) > 0
        j = j + 1
        DadosF(j,1) = Dados(i,6)   
        DadosF(j,2) = Dados(i,7)
        DadosF(j,3) = Dados(i,8)
    end
end
Tm_low = mean(DadosF(:,1))                      // Temperatura média nível inferior [oC]
Tm_up = mean(DadosF(:,2))                       // Temperatura média nível superior [oC]
Um = mean(DadosF(:,3))                         // Velocidade média do escoamento [m/s] 

// Eddy Covariance
for k = 1:j
    Tf(k) = Tm_up-DadosF(k,2)                   // Flutuação da temperatura [oC]
    Wf(k) = Um-DadosF(k,3)                      // Flutuação da velocidade vertical [m/s]
    CovTW(k) = Tf(k)*Wf(k)                      // Covariância TW
end

H_eddy = Ro*Cp*mean(CovTW) 

//Convecção Livre (Campbell e Norman 1997)
wc = 0.81*Lc                                    // Dimensão característica             
rh_cl_cn = 1/(0.05*((Tm_low-Tm_up)/wc)^(1/4))   // Resistência ao transporte de calor 
                                                // convecção livre [m²s/mol]

H_cl_cn = Cpm*(Tm_up-Tm_low)/rh_cl_cn           // Fluxo de calor sensível [W/m²]

//Convecção Forçada (Campbell e Norman 1997)

rh_cf_cn = 7.4*sqrt(wc/Um)                      // Resistência ao transporte de calor 
                                                // convecção forçada [s/m]
H_cf_cn = Cpm*(Tm_up-Tm_low)/rh_cf_cn           // Fluxo de calor sensível [W/m²]

//Convecção Forçada (Monteith and Unsworth 1990)
Re =(Um*Lc)/(Vcf)                                // Número de Reynolds

if Re > 2E4 then    // Fluxo turbulento
    Nu = 0.032*Re^0.8           
    else            // Fluxo laminar
    Nu = 0.60*Re^0.5
end

rh_cf_mu = (Ro*Cp/Dh)*(Lc/Nu)

H_cf_mu = Ro*Cp*(Tm_up-Tm_low)/rh_cf_mu
