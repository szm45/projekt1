## Opis
Program umożliwia konwersję współrzędnych geodezyjnych między układami XYZ, BLH układami PL-1992 i PL-2000. Obsługuje on elipsoidy: WGS84, GRS80 oraz elipsoidę Krasowskiego.
## Wymagania
 Python 3.11,  Biblioteki: argparse, numpy
## System operacyjny
Program został napisany dla systemu operacyjnego Windows 11.
## Instrukcja obsługi
Program umożliwia przekształcenie współrzędnych XYZ na współrzędne BLH, PL-1992 oraz PL-2000 dla wcześniej podanych elipsoid. Dane wejściowe i wyjściowe są obsługiwane w formacie liczby zmiennoprzecinkowej (float). Istnieje możliwość wprowadzenia współrzędnych ręcznie poprzez argument --input cmd.
Aby uzyskać przeliczone współrzędne dla wybranego układu, należy otworzyć wiersz poleceń (cmd) i nawigować do folderu zawierającego plik programu (na przykład: C:\Users\Pulpit\infa\projekt_1). Następnie należy wpisać "python", nazwę pliku programu oraz jedną z wskazanych elipsoid: GRS80, WGS84, Krasowski. Po tej operacji można wprowadzać współrzędne XYZ (wartości współrzędnych wyrażone w metrach) w tej samej linii komend.

*Przykład:* \
*python skrypt.py -m GRS80 -x 57392 -y 9387 -z 4567

*Przykład otrzymanych wyników:* \
Elipsoida: GRS80
Wyniki_hirvonen; fi =  15°06′11.20752″, lam =   9°17′20.44608″, ha = -6319351.650[m]
Wyniki_z_transformacji_1992_oraz_2000; X1992 =  '-' [m], Y1992 =  '-' [m], X2000 =  '-' [m], Y2000 =  '-' [m]
Niewłaściwe położenie

## Znane błędy i nietypowe zachowania
- Program zwraca błąd w przypadku podania niepoprawnego modelu elipsoidy lub systemu współrzędnych.
- Program zwraca błąd dla transformacji XYZ -> BLH  w przypadku podania współrzędnych X=0 Y=0, dla których nie jest możliwe jednoznaczne określenie współrzędnych w układzie BLH. 
