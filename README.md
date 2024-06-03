## Opis
Program umożliwia konwersję współrzędnych geodezyjnych między układami XYZ, FLH układami PL-1992 i PL-2000. Obsługuje on elipsoidy: WGS84, GRS80 oraz elipsoidę Krasowskiego.
## Wymagania
 Python 3.11,  Biblioteki: argparse, numpy
## System operacyjny
Program został napisany dla systemu operacyjnego Windows 11.
## Instrukcja obsługi
Program umożliwia przekształcenie współrzędnych XYZ na współrzędne FLH,FLH na współrzędnę układu PL-1992 oraz PL-2000 dla wcześniej podanych elipsoid. Dane wejściowe i wyjściowe są obsługiwane w formacie liczby zmiennoprzecinkowej (float).

Punktem będącym początkiem układu NEU jest punkt punkt a. Czyli punkt który poprzedza punkt właściwy w pliku tekstowym. A w przypadku pojedyńczego użycia należy podać najpierw środek układu, a następnie obliczany punkt.

Program pozwala na wczytanie współrzędnych z pliku tekstowego (wsp_inp.txt), gdzie dane wejściowe to współrzędne X, Y, Z.Są one oddzielone przecinkiem, a przykałdowy plik znajduje się na githubuie. Przetwarzanie tych danych generuje cztery pliki wynikowe zawierające odpowiednio współrzędne w układach NEU, FLH, X1992Y1992 oraz X2000Y2000. Użytkownik ma możliwość wyboru interesującego go rezultatu i otwarcia odpowiedniego pliku tekstowego.

Program pozwala na wczytanie współrzędnych z pliku tekstowego za pomocą wiersza poleceń. Należy wtedy podać elipsoid, plik z danymi, jednostki pliku wynikowego, nazwy plików wynikowych.

Istnieje możliwość wprowadzenia współrzędnych ręcznie poprzez argument --input cmd. Aby uzyskać przeliczone współrzędne dla wybranego układu, należy otworzyć wiersz poleceń (cmd) i nawigować do folderu zawierającego plik programu (na przykład: C:\Users\Pulpit\infa\projekt_1). Następnie należy wpisać "python", nazwę pliku programu oraz jedną z wskazanych elipsoid: GRS80, WGS84, Krasowski. Po tej operacji można wprowadzać współrzędne XYZ (wartości współrzędnych wyrażone w metrach) w tej samej linii komend.

*Przykład:* \
python skrypt.py -m GRS80 -t input_file.txt -d dms -flh flh123.txt -x92y92 x92y92123.txt -x20y20 x20y20123.txt -neu neu123.txt

Jeśli chcemy przetransformować współrzędne zawarte w pliku tekstowym, wówczas musimy podać elipsoidę (-m GRS80), jego nazwę (-t input_file.txt), jednostkę w pliku wynikowym (-d dms) oraz nazwę plików wyjściowych (-flh flh123.txt -x92y92 x92y92123.txt -x20y20 x20y20123.txt -neu neu123.txt). Ważne jest aby dane w pliku były rozdzielone przecinkiem, a także aby nie zawierał on spacji a współrzędne każdego punktu zaczynały się od nowego wiersza. Separatorem rozwinięcia dziesiętnego liczby powinna być kropka. Należy pamiętać, że plik powinien znajdować się w tym samym folderze roboczym co nasz program.
*Przykład2:* \
*python skrypt.py -m GRS80 -x 57392 -y 9387 -z 4567

*Przykład otrzymanych wyników:* \
Elipsoida: GRS80
Wyniki_hirvonen; fi =  15°06′11.20752″, lam =   9°17′20.44608″, ha = -6319351.650[m]
Wyniki_z_transformacji_1992_oraz_2000; X1992 =  '-' [m], Y1992 =  '-' [m], X2000 =  '-' [m], Y2000 =  '-' [m]
Niewłaściwe położenie

## Znane błędy i nietypowe zachowania
- Program zwraca błąd w przypadku podania niepoprawnego modelu elipsoidy lub systemu współrzędnych.
- Program zwraca błąd dla transformacji XYZ -> BLH  w przypadku podania współrzędnych X=0 Y=0, dla których nie jest możliwe jednoznaczne określenie współrzędnych w układzie BLH.
- Transformacja FLH -> PL-1992 oraz FLH -> PL-2000 dla elipsoidy krasowskiego wzraca błędne wyniki, dlatego nie można ich używać
