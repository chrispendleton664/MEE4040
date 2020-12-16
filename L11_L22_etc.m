L11*u0 + L12*v0 + L13*w0 + N1*w0 + mdot*u02 - qx = 0
L21u0 + L22*v0 + L23*w0 + N2*w0 + mdot*v02 - qy = 0
L31u0 + L32*v0 + L33*w0 + N2*(u0, v0, w0?) + mdot*w02 - qz = 0

% L1j
% X
x = sym("x");
A11x = Aijsum(1,1);
A1101x = diff(A11x,x);
A1102x = diff(A1101x,x);
A13x = Aijsum(1,3);
A1301x = diff(A13x,x);
A1302x = diff(A1301x,x);
A33x = Aijsum(3,3);
A3301x = diff(A33x,x);
A3302x = diff(A3301x,x);

B11x = Bijsum(1,1);
B1101x = diff(B11x,x);
B1102x = diff(B1101x,x);
B1103x = diff(B1102x,x);
B13x = Bijsum(1,3);
B1301x = diff(B13x,x);
B1302x = diff(B1301x,x);
B1303x = diff(B1302x,x);

D11x = Dijsum(1,1);
D1101x = diff(D11x,x);
D1102x = diff(D1101x,x);
D1103x = diff(D1102x,x);
D1104x = diff(D1103x,x);

% Y
y = sym("y");
A11y = Aijsum(1,1);
A1101y = diff(A11y,y);
A1102y = diff(A1101y,y);
A22y = Aijsum(2,2);
A2201y = diff(A22y,y);
A2202y = diff(A2201y,y);
A23y = Aijsum(2,3);
A2301y = diff(A23y,y);
A2302y = diff(A2301y,y);
A33y = Aijsum(3,3);
A3301y = diff(A33y,y);
A3302y = diff(A3301y,y);

B22y = Bijsum(1,1);
B2201y = diff(B22y,y);
B2202y = diff(B2201y,y);
B2203y = diff(B2202y,y);
B23y = Bijsum(1,1);
B2301y = diff(B23y,y);
B2302y = diff(B2301y,y);
B2303y = diff(B2302y,y);

D22y = Dijsum(2,2);
D2201y = diff(D23y,y);
D2202y = diff(D2301y,y);
D2203y = diff(D2302y,y);
D2204y = diff(D2303y,y);

% XY
A12xy = Aijsum(1,2);
A1201xy = diff(A12xy,x);
A1202xy = diff(A1201xy,y);
A13xy = Aijsum(1,3);
A1301xy = diff(A13xy,x);
A1302xy = diff(A1301xy,y);
A23xy = Aijsum(2,3);
A2301xy = diff(A23xy,x);
A2302xy = diff(A2301xy,y);
A33xy = Aijsum(3,3);
A3301xy = diff(A33xy,x);
A3302xy = diff(A3301xy,y);

B12xxy = Bijsum(1,2);
B1201xxy = diff(B12xxy,x);
B1202xxy = diff(B1201xxy,x);
B1203xxy = diff(B1202xxy,y);
B13xxy = Bijsum(1,3);
B1301xxy = diff(B13xxy,x);
B1302xxy = diff(B1301xxy,x);
B1303xxy = diff(B1302xxy,y);
B33xxy = Bijsum(3,3);
B3301xxy = diff(B33xxy,x);
B3302xxy = diff(B3301xxy,x);
B3303xxy = diff(B3302xxy,y);
B12xyy = Bijsum(1,2);
B1201xyy = diff(B12xyy,x);
B1202xyy = diff(B1201xyy,y);
B1203xyy = diff(B1202xyy,y);
B23xyy = Bijsum(2,3);
B2301xyy = diff(B23xyy,x);
B2302xyy = diff(B2301xyy,y);
B2303xyy = diff(B2302xyy,y);
B33xyy = Bijsum(3,3);
B3301xyy = diff(B33xyy,x);
B3302xyy = diff(B3301xyy,y);
B3303xyy = diff(B3302xyy,y);

D11xxxy = Dijsum(1,1);
D1101xxxy = diff(D11xxxy,x);
D1102xxxy = diff(D1101xxxy,x);
D1103xxxy = diff(D1102xxxy,x);
D1104xxxy = diff(D1103xxxy,y);
D12xxyy = Dijsum(1,2);
D1201xxyy = diff(D12xxyy,x);
D1202xxyy = diff(D1201xxyy,x);
D1203xxyy = diff(D1202xxyy,y);
D1204xxyy = diff(D1203xxyy,y);
D33xxyy = Dijsum(1,2);
D3301xxyy = diff(D33xxyy,x);
D3302xxyy = diff(D3301xxyy,x);
D3303xxyy = diff(D3302xxyy,y);
D3304xxyy = diff(D3303xxyy,y);
D23xyyy = Dijsum(2,3);
D2301xyyy = diff(D23xyyy,x);
D2302xyyy = diff(D2301xyyy,y);
D2303xyyy = diff(D2302xyyy,y);
D2304xyyy = diff(D2303xyyy,y);


