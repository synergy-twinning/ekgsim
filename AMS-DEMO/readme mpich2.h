1. inštalira se mpich2 (package manager v Ubuntu)
2. nastavi se mpd
	mpd.hosts v domačem direktoriju uporabnika naj vsebuje seznam nodov (ločenih z novo vrstico) in število procesorjev na njih (npr: localhost:8)
	.mpd.conf v domačem direktoriju uporabnika vsebuje definicijo šifre (ne vem zakaj se uporablja) (secretword=<izbrana šifra>)
	požene se mpdboot , ki se mu poda število hostov (mpdboot 1)
3. računalniki morajo imeti nastavljen ssh strežnik (ubuntu ga po defaultu nima - instalira se ssh paket)
	na enem računalniku se zgradi nezaščiten dsa ključ, katerega pub vrednost se zapiše v ~/.ssh/known_hosts vseh sodelujočih računalnikov

