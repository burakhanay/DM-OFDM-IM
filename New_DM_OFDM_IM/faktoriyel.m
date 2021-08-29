function sonuc= faktoriyel(sayi)
sonuc = sayi;
if  sayi == 0
    sonuc = 1;
else
    sonuc = sonuc *  faktoriyel(sayi-1);
end