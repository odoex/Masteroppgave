# Masteroppgave

Den stabile mesh refinement metoden fungerer nå for Advection equation. (Tror jeg)
Den ligger under Advection_equation\stableMethod. Tilhørende filer er: 
-meshRefinementAdv2.m
-rhsFineGridAdv.m
-projectFineBoundary.m

Tilhørende filer som også brukes i vanlig mesh refinement advection er:
-createSolutionVector.m
-finiteVolumeAdv.m
-RK_4Adv.m
-exactSolAdv.m
-f_adv.m
-g_adv.m
-rhsAdv.m
-boundary.m

De to førstenevnte, finiteVolumeAdv.m og RK_4Adv.m,
må modereres for å bruke den stabile metoden. For finiteVolumeAdv.m 
byttes den vanlige projectFlux() ut med projectFineBoundary() (rart navn, 
skal byttes). For RK_4.m kommenteres den vanlige rk-metoden ut og den 
nedenforstående metoden brukes, hvor det sjekkes om gridet er et undergrid
hvor det i så fall kalles på funksjonen rhsFineGridAdv() som håndterer 
overgangen mellom fint og grovt grid på en annen måte enn som en rand. To 
av de tre testene i mappen tests fungerer ikke selv om det er riktig. Må
finne ut av hvorfor.

Euler equations
Euler ligningene sin mesh refinement metode fungerer så vidt jeg vet, men 
gir ikke ønskede resultater ettersom den har  mye lavere konvergens. Den 
befinner seg i mappen Euler_equations. 
