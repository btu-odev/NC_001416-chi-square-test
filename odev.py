import re
import math
from Bio import SeqIO
from scipy.integrate import quad
from terminaltables import AsciiTable

pathToFile = "NC_001416.fasta" # Canlımızın sekanslarının bulunduğu dosya

fasta_sequences = SeqIO.parse(open(pathToFile), 'fasta') # Dosyamızdaki tüm fasta sekanslarını getiren Iterator
sequence_obj = next(fasta_sequences) # Iteratorü bir adım ileriye götürerek canlımızın nesnesine ulaşıyoruz

id, sequence = sequence_obj.id, str(sequence_obj.seq) # Nesnemizden id ve sekansı alıyoruz

HpaII = 'CCGG' # Kullanacağımız enzimin parçaladığı sekans
h = [m.start() for m in re.finditer(HpaII, sequence)] # HpaII içeren stringlerin başlangıç indislerini tutuyor

fragments = ["" for x in range(len(h))] # Fragmanlar için boz bir liste oluşturuyoruz

for i in range(len(h) - 1):
    fragments[i] = sequence[h[i]:h[i+1] - 1] # Enzimin böleceği noktalardan fragmanları elde ediyoruz
fragments[-1] = sequence[h[-1]:] # Son fragmanı elde ediyoruz

p_C2 = (len([i.start() for i in re.finditer('C', sequence)]) / (len(sequence)))**2 # CC dimerinin bulunma sıklığı
p_G2 = (len([i.start() for i in re.finditer('G', sequence)]) / (len(sequence)))**2 # GG dimerinin bulunma sıklığı
p = p_C2*p_G2 # CCGG sekansının p-value değeri

poisson = lambda x: p*math.exp(-p*x) # Poisson dağılımını hesaplayan lambda fonksiyonumuz
bins = list(range(0, 601, 100)) + [len(sequence)] # Aralıklarımızın tutulduğu dizi

expecteds = [] # Beklenen değerlerimizin tutulduğu dizi
for i in range(len(bins) - 1):
    res, err = quad(poisson, bins[i], bins[i+1]) # Belirli integral hesaplanıyor
    res = round(res, 3)
    expecteds.append(len(fragments)*res) # Hesaplanan olasılık dizinin sonuna ekleniyor

observations = [0 for i in range(len(bins) - 1)] # Her aralık için bir hücre oluşturuluyor
for i in fragments:
    for j in range(len(bins) - 1):
        if len(i) >= bins[j] and len(i) < bins[j+1]: # Fragman uzunluğu aralığa dahilse ilgili kovanın değerini artırıyoruz
            observations[j] += 1
            break

chi_square_test = [] # Ki-kare testlerimizin tutulduğu dizi
for expected, observation in zip(expecteds, observations):
    chi_square_test.append((observation - expected) ** 2 / expected)

# Tablo formatına uyulabilmesi için diziler düzenleniyor

bins = ['Aralıklar'] + [str(bins[i]) + ' - ' + str(bins[i+1]) for i in range(len(bins) - 1  )] + ['------'] + ['Toplam']
observations = ['Gözlemler'] + observations + ['------'] + [round(sum(observations), 3)]
expecteds = ['Beklenen'] + [round(i, 3) for i in expecteds] + ['------'] + [round(sum(expecteds), 3)]
chi_square_test = ['Ki-kare Test'] + [round(i, 3) for i in chi_square_test] + ['------'] + [round(sum(chi_square_test), 3)]

# Oluşturulan diziler ile terminal üzerinde tablo oluşturuluyor
table_data = [[a,b,c,d] for a,b,c,d in zip(bins, observations, expecteds, chi_square_test)]
table = AsciiTable(table_data)
print(table.table)
