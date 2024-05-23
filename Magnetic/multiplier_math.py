# Define the inverse cube law function
def Magnet(x, A, B):
    y = B/((x+A)**3)
    return y

def ArrayForce(d, n, A, B):
	F = Magnet(d, A, B) # [lb]
	Ftot = F*n 
	return Ftot

A = 0.2870475701906723 # [in]
B = 0.7334768699294826 # [lb*in^3] Treating this as inverse cube law | F = B/(x+A)^3, will consult with Zack for better modeling
n = 30
d = 0.4 # [in]

Ftot = ArrayForce(d, n, A, B)
print(Ftot) # Compared with real value of 135 lb | K*J gave 64.8 lb at 0.4in
print(135)

print(135/Ftot)

print(A*2.54)
print(B*2.54**3/0.225)