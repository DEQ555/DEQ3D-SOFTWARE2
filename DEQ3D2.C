#include <stdlib.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////
// memset
static void _Memset( void* _x, char _y, unsigned _z ){
	char *_X_=(char*)_x;
	for (unsigned int i = 0; i < _z; ++i)_X_[i] = _y;
}

////////////////////////////////////////////////////////////////////////////
// Save file .bmp
void __fastcall SaveBMP( const char* filename, unsigned char *data, const unsigned int ancho, const unsigned int alto, unsigned char bpp ){
	int Y;
	unsigned int X,STEP;
	unsigned char BPP, BYTES_BAZURAS[3];
	if ( bpp<24 )bpp*=8;
	if ( !data ) return;
	const unsigned char padding = ( ( 4 - ( ancho * (bpp/8) ) % 4) % 4 );
	unsigned char header[14], informacionheader[ 40 ];
	_Memset( header, 0, ( sizeof(unsigned char) * 14 ) );
	_Memset( informacionheader, 0, ( sizeof(unsigned char) * 40 ) );
	const unsigned int filesize = 14 + 40 + ancho * alto * (bpp/8) + padding * ancho;
	header[0] = 'B',header[1] = 'M';
	header[2] = filesize, header[3] = filesize >> 8, header[4] = filesize >> 16, header[5] = filesize >> 24;
	header[10] = 14 + 40;
	informacionheader[0] = 40;
	informacionheader[4] = ancho, informacionheader[5] = ancho >> 8, informacionheader[6] = ancho >> 16, informacionheader[7]  = ancho >> 24;
	informacionheader[8] = alto,  informacionheader[9] = alto >> 8, informacionheader[10] = alto >> 16,  informacionheader[11] = alto >> 24;
	informacionheader[12] = 1, informacionheader[14] = bpp;
	FILE *f = fopen( filename, "wb" );
	fwrite( header, 14, 1, f );
	fwrite( informacionheader, 40, 1, f );
	BPP = (bpp/8);
	for ( Y = alto-1; Y >=0 ; --Y){
		for ( X = 0; X < ancho; ++X) STEP = ( Y * ancho + X ) * BPP, fwrite( &data[ STEP ], BPP, 1, f );
		fwrite( BYTES_BAZURAS, padding, 1, f );
	}
	fclose(f),f=NULL;
	return;
}




////////////////////////////////////////////////
// 
//
typedef char _s8;
typedef short _s16;
typedef int _s32;
typedef unsigned char _u8;
typedef unsigned short _u16;
typedef unsigned int _u32;

struct GUB_FVM3D_HEAD{
	unsigned char     Modo; // En cada bit en binario representara algun modo
	unsigned int Num_Tam_V; // Numero de multiplicado Vertices
	unsigned int Num_Tam_T; // Numero de multiplicado Vertices Textura
	unsigned int     Num_V; // Numero de Vertices
	unsigned int     Num_T; // Numero de Vertices Textura
	unsigned int     Num_F; // Numero de Faces
	unsigned int   Num_TEX; // Numero de Texturas
};

struct GUB_FVM3D_OBJECT{
	// Estructura de Cabesera FVM
	struct GUB_FVM3D_HEAD FVM_HEAD;
	unsigned int *TEX; // Buffer Texturas
	int           *V_; // Buffer Vertices
	unsigned int  *T_; // Buffer Vertices Textura
	unsigned int  *F_; // Buffer Faces
};


void _GUB_NULL_FVM_( struct GUB_FVM3D_OBJECT *X ){
	X->FVM_HEAD.Modo = 0;
	X->FVM_HEAD.Num_Tam_V = 0;
	X->FVM_HEAD.Num_Tam_T = 0;
	X->FVM_HEAD.Num_V = 0;
	X->FVM_HEAD.Num_T = 0;
	X->FVM_HEAD.Num_F = 0;
	X->FVM_HEAD.Num_TEX = 0;
	X->V_ = NULL;
	X->T_ = NULL;
	X->F_ = NULL;
	X->TEX = NULL;
}

void _GUB_FREE_FVM_( struct GUB_FVM3D_OBJECT *X ){
	free(X->V_),  X->V_  = NULL;
	free(X->T_),  X->T_  = NULL;
	free(X->F_),  X->F_  = NULL;
	free(X->TEX), X->TEX = NULL;
	X->FVM_HEAD.Modo = X->FVM_HEAD.Num_Tam_V = X->FVM_HEAD.Num_Tam_T = X->FVM_HEAD.Num_V = X->FVM_HEAD.Num_T = X->FVM_HEAD.Num_F = X->FVM_HEAD.Num_TEX = 0;
}

int _GUB_LOAD_FVM_( struct GUB_FVM3D_OBJECT *X, const char Archivo[] ){
	FILE *f;
	f = fopen( Archivo, "rb" );
	fread( &X->FVM_HEAD, sizeof(struct GUB_FVM3D_HEAD), 1, f );

	X->TEX = (unsigned int*)malloc( (sizeof(unsigned int) * X->FVM_HEAD.Num_Tam_T * X->FVM_HEAD.Num_Tam_T) * X->FVM_HEAD.Num_TEX );
	X->V_ = (int*)malloc( sizeof(int) * X->FVM_HEAD.Num_V * 3);
	X->T_ = (unsigned int*)malloc( sizeof(unsigned int) * X->FVM_HEAD.Num_T * 2 );
	X->F_ = (unsigned int*)malloc( sizeof(unsigned int) * X->FVM_HEAD.Num_F * 8 );

	fread( X->TEX, (sizeof(unsigned int) * X->FVM_HEAD.Num_Tam_T * X->FVM_HEAD.Num_Tam_T) * X->FVM_HEAD.Num_TEX, 1, f );
	fread( X->V_, sizeof(int) * X->FVM_HEAD.Num_V * 3, 1, f );
	fread( X->T_, sizeof(unsigned int) * X->FVM_HEAD.Num_T * 2, 1, f );
	fread( X->F_, sizeof(unsigned int) * X->FVM_HEAD.Num_F * 8, 1, f );
	fclose(f);
	f=NULL;
	return 0;
}

////////////////////////////////////////////////
// Grill Draw
//
typedef struct{
	float area;
	float vx0,vy0,vz0;
	float vx1,vy1,vz1;
	float vx2,vy2,vz2;
}GRILL;

////////////////////////////////////////////////
// Create Screen Buffer[CANVAS]
//
int* Vbuff = NULL;
int* Zbuff = NULL;
GRILL* Gbuff = NULL;
unsigned int _size = 0;
unsigned int _w = 0;
unsigned int _h = 0;
unsigned int _wm = 0;
unsigned int _hm = 0;

struct GUB_FVM3D_OBJECT _Model;

void CreateCanvas(
	const int w,
	const int h
){
	_w = w;
	_h = h;
	_wm = w>>1;
	_hm = h>>1;
	_size = w * h;
	Vbuff = (int*)malloc( sizeof(int) * _size );
	Zbuff = (int*)malloc( sizeof(int) * _size );
	Gbuff = (GRILL*)malloc( sizeof(GRILL) * _size );
	for (unsigned int i = 0; i < _size; ++i){
		Vbuff[i]=0;
		Zbuff[i]=0;
		Gbuff[i].area=0;
	}
}

void FreeCanvas(void){
	_wm=_hm=_w=_h=_size=0;
	free(Vbuff);
	free(Zbuff);
	free(Gbuff);
	Vbuff=NULL;
	Zbuff=NULL;
	Gbuff=NULL;
}

////////////////////////////////////////////////
// Draw Screen/Canvas
//
  void __fastcall FlipBuff(void){
	SaveBMP("VideoRaster.bmp",(unsigned char*)Vbuff,_w,_h,32);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/*          RASTER          */

// Función para dibujar un triángulo en modo TriangleWire
/*void DrawTriangleWire(int x0, int y0, int x1, int y1, int x2, int y2, int color) {
    // Dibuja las tres líneas del triángulo que conectan sus vértices
    int dx01 = abs(x0 - x1), dy01 = abs(y0 - y1);
    int dx12 = abs(x1 - x2), dy12 = abs(y1 - y2);
    int dx20 = abs(x2 - x0), dy20 = abs(y2 - y0);

    int longest = dx01;
    if (dx12 > longest) longest = dx12;
    if (dx20 > longest) longest = dx20;

    for (int i = 0; i <= longest; i++) {
        float t01 = (float)i / longest;
        float t12 = (float)i / longest;
        float t20 = (float)i / longest;

        int x01 = x0 + dx01 * t01, y01 = y0 + dy01 * t01;
        int x12 = x1 + dx12 * t12, y12 = y1 + dy12 * t12;
        int x20 = x2 + dx20 * t20, y20 = y2 + dy20 * t20;

        PixelSet(x01, y01, color);
        PixelSet(x12, y12, color);
        PixelSet(x20, y20, color);
    }
}

void DrawTriangleWireOptimized(int x0, int y0, int x1, int y1, int x2, int y2, int* Vbuff, int Ancho, int Alto, int color) {
    // Función para establecer un píxel en el búfer RGBA
    // Supongamos que cada píxel se almacena como un entero RGBA en Vbuff.
    // El formato RGBA se interpreta como: R en los primeros 8 bits, G en los siguientes 8, B en los siguientes 8 y A en los últimos 8.

    // Clampea las coordenadas dentro de los límites de la pantalla
    x0 = (x0 < 0) ? 0 : (x0 >= Ancho) ? Ancho - 1 : x0;
    y0 = (y0 < 0) ? 0 : (y0 >= Alto) ? Alto - 1 : y0;
    x1 = (x1 < 0) ? 0 : (x1 >= Ancho) ? Ancho - 1 : x1;
    y1 = (y1 < 0) ? 0 : (y1 >= Alto) ? Alto - 1 : y1;
    x2 = (x2 < 0) ? 0 : (x2 >= Ancho) ? Ancho - 1 : x2;
    y2 = (y2 < 0) ? 0 : (y2 >= Alto) ? Alto - 1 : y2;

    // Encuentra las diferencias y direcciones para cada borde del triángulo
    int dx01 = abs(x0 - x1), dy01 = abs(y0 - y1);
    int dx12 = abs(x1 - x2), dy12 = abs(y1 - y2);
    int dx20 = abs(x2 - x0), dy20 = abs(y2 - y0);

    int sx01 = (x0 < x1) ? 1 : -1, sy01 = (y0 < y1) ? 1 : -1;
    int sx12 = (x1 < x2) ? 1 : -1, sy12 = (y1 < y2) ? 1 : -1;
    int sx20 = (x2 < x0) ? 1 : -1, sy20 = (y2 < y0) ? 1 : -1;

    // Calcula el punto de inicio y las direcciones de los bordes
    int x = x0, y = y0;
    int x1a = x1, y1a = y1;
    int x2a = x2, y2a = y2;
    int dx1, dy1, dx2, dy2;

    if (dy01 > dx01) {
        dx1 = dy01;
        dy1 = dx01;
        x1a = x0 + sx01;
    } else {
        dx1 = dx01;
        dy1 = dy01;
        y1a = y0 + sy01;
    }

    if (dy12 > dx12) {
        dx2 = dy12;
        dy2 = dx12;
        x2a = x1 + sx12;
    } else {
        dx2 = dx12;
        dy2 = dy12;
        y2a = y1 + sy12;
    }

    // Dibuja los bordes del triángulo
    while ((x != x1a || y != y1a) && (x != x2a || y != y2a)) {
        // Dibuja el píxel actual
        if (x >= 0 && x < Ancho && y >= 0 && y < Alto) {
            int index = y * Ancho + x;
            Vbuff[index] = color; // Establece el color en el búfer RGBA
        }

        int e1 = 2 * err1;
        int e2 = 2 * err2;

        if (e1 >= dy1) {
            err1 += dy1;
            x += sx01;
        }

        if (e2 >= dy2) {
            err2 += dy2;
            x += sx12;
        }

        if (e1 <= dx1) {
            err1 += dx1;
            y += sy01;
        }

        if (e2 <= dx2) {
            err2 += dx2;
            y += sy12;
        }
    }
}

*/

#define _SWAP(A,B){ A^=B,B^=A,A^=B; }
typedef struct{int x,y;}PUNTOS;

void LineW(int x0, int x1, int y0, int color){
	if(x0>x1)_SWAP(x0,x1);
	int Len = (x1-x0);
	while(Len--){
		Vbuff[ y0 * _w + x0++ ] = color;
	}
}

void polygonfill(const unsigned int cantidad, PUNTOS p[], int color) {
    int xmin = p[0].x;
    int xmax = p[0].x;
    int ymin = p[0].y;
    int ymax = p[0].y;

    unsigned int i;
  
    for (i = 1; i < cantidad; i++) {
        if (xmin > p[i].x) xmin = p[i].x;
        if (xmax < p[i].x) xmax = p[i].x;
        if (ymin > p[i].y) ymin = p[i].y;
        if (ymax < p[i].y) ymax = p[i].y;
    }
  
    double s = ymin + 0.01;
  
    while (s <= ymax) {
        int inter[cantidad], c = 0;
        int x1, x2, y1, y2, x, *Ti;
      
        for (i = 0; i < cantidad; i++) {
            x1 = p[i].x;
            y1 = p[i].y;
            x2 = p[(i + 1) % cantidad].x;
            y2 = p[(i + 1) % cantidad].y;
          
            if (y2 < y1) {
                _SWAP(x1, x2);
                _SWAP(y1, y2);
            }
          
            if (s <= y2 && s >= y1) {
                if ((y1 - y2) == 0)
                    inter[c++] = x1;
                else {
                    x = ((x2 - x1) * (s - y1)) / (y2 - y1);
                    inter[c++] = x + x1;
                }
            }
        }
        Ti=inter;
        while(c){
        	LineW(*Ti++,*Ti++,s,color);
        	c-=2;
        }
        s++;
    }
}

#define __AREA(X0,Y0,X1,Y1,X2,Y2) (float)(((X2) - (X0)) * ((Y1) - (Y0)) - ((Y2) - (Y0)) * ((X1) - (X0)))
#define __near 300
#define __fac 900 * (__near)
void __fastcall DrawModel(float Cam[3]){
	unsigned int i;
	int _x,_y,_z,FOV;
	PUNTOS Points[4];
	int _V[_Model.FVM_HEAD.Num_V<<2];
	int* Pv = _Model.V_;
	unsigned int* Pf = _Model.F_;
	unsigned int index = 0;
	for (i=0; i<_Model.FVM_HEAD.Num_V; ++i){
		_x = *Pv++,
		_y = *Pv++,
		_z = *Pv++;

		_x += Cam[0];
		_y += Cam[1];
		_z += Cam[2];

		FOV = __fac / (_z+(1>>_z));
		_V[index+0] = ( ( ( _x * FOV) ) >> 10 ) + _wm;
		_V[index+1] = ( ( ( _y * FOV) ) >> 10 ) + _hm;
		_V[index+2] = FOV;
		_V[index+3] = _z;

		Vbuff[ _V[index+1]*_w+_V[index+0] ] = 0xffffffff;

		index += 4;
	}
	for (i=0; i<_Model.FVM_HEAD.Num_F; ++i){
		int Av = (*Pf++)<<2;
		int Bv = (*Pf++)<<2;
		int Cv = (*Pf++)<<2;

		int At = (*Pf++)<<1;
		int Bt = (*Pf++)<<1;
		int Ct = (*Pf++)<<1;
		Pf++;
		Pf++;

		Points[0].x=_V[Av++];
		Points[0].y=_V[Av];
		Points[1].x=_V[Bv++];
		Points[1].y=_V[Bv];
		Points[2].x=_V[Cv++];
		Points[2].y=_V[Cv];

		float A0 = __AREA(
			Points[0].x,
			Points[0].y,
			
			Points[1].x,
			Points[1].y,
			
			Points[2].x,
			Points[2].y
		);

		if(A0>=0)continue;

		polygonfill(3,Points,0xffff0000);
	}
}
void __fastcall SetShader(void){
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


int main(int argc, char const *argv[]){

	_GUB_NULL_FVM_(&_Model);
	// _GUB_LOAD_FVM_(&_Model,"extra.fvm");
	_GUB_LOAD_FVM_(&_Model,"soldier.fvm");

	CreateCanvas(640,480);

	float Camera[3] = {0,0,4000};

	DrawModel(Camera);
	SetShader();

	FlipBuff();
	FreeCanvas();

	_GUB_FREE_FVM_(&_Model);

	return 0;
}