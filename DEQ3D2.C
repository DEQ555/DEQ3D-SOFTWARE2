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
void SaveBMP( const char* filename, unsigned char *data, const unsigned int ancho, const unsigned int alto, unsigned char bpp ){
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
	int act;
	float area;
	unsigned int Pmodel;
	float vx0,vy0,vz0;
	float vx1,vy1,vz1;
	float vx2,vy2,vz2;
	float tx0,ty0;
	float tx1,ty1;
	float tx2,ty2;
}GRILL;
#include <intrin.h>

////////////////////////////////////////////////
// Create Screen Buffer[CANVAS]
//
int* Vbuff = NULL;
int* Zbuff = NULL;
GRILL* Gbuff = NULL;
int _fills[2] = {0,-99999};
unsigned int _size = 0;
unsigned int _w = 0;
unsigned int _h = 0;
unsigned int _wm = 0;
unsigned int _hm = 0;

struct GUB_FVM3D_OBJECT _Model[5];

void CreateCanvas(
	const int w,
	const int h
){
	_w = w;
	_h = h;
	_wm = w>>1;
	_hm = h>>1;
	_size = w * h;
	Vbuff = (int*)_aligned_malloc( sizeof(int) * _size, sizeof(unsigned int) * 8 );
	Zbuff = (int*)_aligned_malloc( sizeof(int) * _size, sizeof(unsigned int) * 8 );
	Gbuff = (GRILL*)_aligned_malloc( sizeof(GRILL) * _size, sizeof(unsigned int) * 8 );
	for (unsigned int i = 0; i < _size; ++i){
		Vbuff[i]=_fills[0];
		Zbuff[i]=_fills[1];
		Gbuff[i].act=0;
	}
}

void FreeCanvas(void){
	_wm=_hm=_w=_h=_size=0;
	_aligned_free(Vbuff);
	_aligned_free(Zbuff);
	_aligned_free(Gbuff);
	Vbuff=NULL;
	Zbuff=NULL;
	Gbuff=NULL;
}

void SetFillCanvas( int c ){
	_fills[0]=c;
}

void SetFillDepth( int c ){
	_fills[1]=c;
}

////////////////////////////////////////////////
// Draw Screen/Canvas
//
  void FlipBuff(void){
	SaveBMP("VideoRaster.bmp",(unsigned char*)Vbuff,_w,_h,32);
}




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

#define MAXSIN 255
const unsigned char sinTab[91] = {
0,4,8,13,17,22,26,31,35,39,44,48,53,57,61,65,70,74,78,83,87,91,95,99,103,107,111,115,119,123,
127,131,135,138,142,146,149,153,156,160,163,167,170,173,177,180,183,186,189,192,195,198,200,203,206,208,211,213,216,218,
220,223,225,227,229,231,232,234,236,238,239,241,242,243,245,246,247,248,249,250,251,251,252,253,253,254,254,254,254,254,
255
};

int fastSin(int i)
{
	while(i<0) i+=360;
	while(i>=360) i-=360;
	if(i<90)  return((sinTab[i])); else
	if(i<180) return((sinTab[180-i])); else
	if(i<270) return(-(sinTab[i-180])); else
	return(-(sinTab[360-i]));
}

int fastCos(int i){
	return fastSin(i+90);
}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/*          RASTER          */

#define _SWAP(A,B){ A^=B,B^=A,A^=B; }
typedef struct{int x,y;}PUNTOS;

#define __AREA(X0,Y0,X1,Y1,X2,Y2) (float)(((X2) - (X0)) * ((Y1) - (Y0)) - ((Y2) - (Y0)) * ((X1) - (X0)))
#define __near 300
#define __fac 900 * (__near)
void DrawModel(const unsigned int Pmodel, float Cam[3],float Rot[3]){
	unsigned int i;
	int _x,_y,_z,FOV;
	GRILL Att;
	int _V[ _Model[Pmodel].FVM_HEAD.Num_V<<2 ];
	int* Pv = _Model[Pmodel].V_;
	unsigned int* Pt = _Model[Pmodel].T_;
	unsigned int* Pf = _Model[Pmodel].F_;
	unsigned int index = 0;
	/*#####################*/
	int cos0,sin0,cos1,sin1;
	cos0 = fastCos(Rot[0]);
	sin0 = fastSin(Rot[0]);
	cos1 = fastCos(Rot[1]);
	sin1 = fastSin(Rot[1]);
	/*#####################*/
	for (i=0; i<_Model[Pmodel].FVM_HEAD.Num_V; ++i){
		int _x_ = *Pv++,
			_y_ = *Pv++,
			_z_ = *Pv++;

		_x = (cos0*_x_ + sin0*_z_)/MAXSIN;
		_y = (cos1*_y_ + (cos0*sin1*_z_-sin0*sin1*_x_)/MAXSIN)/MAXSIN;
		_z = ((cos0*cos1*_z_-sin0*cos1*_x_)/MAXSIN - sin1*_y_)/MAXSIN;

		_x += Cam[0];
		_y += Cam[1];
		_z += Cam[2];

		FOV = __fac / (_z+(1>>_z));
		_V[index+0] = ( ( ( _x * FOV) ) >> 10 ) + _wm;
		_V[index+1] = ( ( ( _y * FOV) ) >> 10 ) + _hm;
		_V[index+2] = FOV;
		_V[index+3] = _z;
		index += 4;
	}
	for (i=0; i<_Model[Pmodel].FVM_HEAD.Num_F; i++){
		unsigned int Av = (*Pf++)*4;
		unsigned int Bv = (*Pf++)*4;
		unsigned int Cv = (*Pf++)*4;
		
		unsigned int At = (*Pf++)*2;
		unsigned int Bt = (*Pf++)*2;
		unsigned int Ct = (*Pf++)*2;
		Pf+=2;

		Att.vx0 = _V[ Av++ ];
		Att.vy0 = _V[ Av++ ];
		Att.vz0 = _V[ Av ];
		Att.vx1 = _V[ Bv++ ];
		Att.vy1 = _V[ Bv++ ];
		Att.vz1 = _V[ Bv ];
		Att.vx2 = _V[ Cv++ ];
		Att.vy2 = _V[ Cv++ ];
		Att.vz2 = _V[ Cv ];

		Att.tx0 = Pt[ At++ ] / Att.vz0;
		Att.ty0 = Pt[ At ] / Att.vz0;
		Att.tx1 = Pt[ Bt++ ] / Att.vz1;
		Att.ty1 = Pt[ Bt ] / Att.vz1;
		Att.tx2 = Pt[ Ct++ ] / Att.vz2;
		Att.ty2 = Pt[ Ct ] / Att.vz2;

		Att.area = __AREA( Att.vx0, Att.vy0, Att.vx1, Att.vy1, Att.vx2, Att.vy2);
		Att.act = 0;
		Att.Pmodel = Pmodel;
		if(
			(unsigned int)Att.vz0 >= __fac ||
			(unsigned int)Att.vz1 >= __fac ||
			(unsigned int)Att.vz2 >= __fac
		)continue;
		if(Att.area>=0)continue;
		{
			int min_y = __min( Att.vy0, __min( Att.vy1, Att.vy2 ) );
			int min_x = __min( Att.vx0, __min( Att.vx1, Att.vx2 ) );
			int max_y = __max( Att.vy0, __max( Att.vy1, Att.vy2 ) );
			int max_x = __max( Att.vx0, __max( Att.vx1, Att.vx2 ) );
			if(min_y<0){ min_y=0; }
			if(max_y>=_h){ max_y=_h-1; }
			if(min_x<0){ min_x=0; }
			if(max_x>=_w){ max_x=_w-1; }
			_z = (( Att.vz0 + Att.vz1 + Att.vz2 ) / 3);
			Att.act = 1;
			for (_y=min_y; _y<max_y; ++_y)
			{
				for (_x=min_x; _x<max_x; ++_x)
				{
					float W0 = __AREA( Att.vx1,Att.vy1,Att.vx2,Att.vy2,_x,_y)/Att.area;
					float W1 = __AREA( Att.vx2,Att.vy2,Att.vx0,Att.vy0,_x,_y)/Att.area;
					float W2 = __AREA( Att.vx0,Att.vy0,Att.vx1,Att.vy1,_x,_y)/Att.area;
					if ( (Zbuff[ _y * _w + _x ] < _z) && ((W0>=0) & (W1>=0) & (W2>=0)) ){
						Zbuff[ _y * _w + _x ] = _z;
						Gbuff[ _y * _w + _x ] = Att;
					}
				}
			}
		}
	}
}

int ___RGB(int r, int g, int b){
	return ( (r&0xff) << 16 | (g&0xff)<<8 | (b&0xff) );
}

int BilinealTex(const unsigned int X, const unsigned int Y, const unsigned int Tam, unsigned int* Tex){
	float totalR = 0, totalG = 0, totalB = 0;
	unsigned int count=0,r, c;
	unsigned char *channel;
	const unsigned int xaxis[] = {Y - 1, Y, Y + 1}, yaxis[] = {X - 1, X, X + 1};
	for ( r = 0; r < 3; ++r){
	    for ( c = 0; c < 3; ++c){
			const int curRow = xaxis[r];
			const int curCol = yaxis[c];
	    	if ( (unsigned int)curRow < Tam && (unsigned int)curCol < Tam ){
	    		channel = (unsigned char*)&Tex[ curRow * Tam + curCol ];
	    		totalB += *channel++;
	    		totalG += *channel++;
	    		totalR += *channel++;
	    		count++;
	    	}
	    }
	}
	totalR /= count, totalG /= count, totalB /= count;
	return ___RGB(totalR,totalG,totalB);
}

#define _Efect_RGB(_X,_Y) ((X+Y-X+Y)-X | (X+Y+X+Y)-Y * 0xff)

void SetShader(void){
	int X,Y,U,V;
	float W0,W1,W2,Z;
	int* VIDEO = Vbuff;
	int* DEPTH = Zbuff;
	GRILL* _X_ = Gbuff;
	for (Y = 0; Y < _h; ++Y)
	{
		for (X = 0; X < _w; ++X)
		{
			*VIDEO = _fills[0] & _Efect_RGB(X,Y);
			*DEPTH = _fills[1];
			W0 = __AREA( _X_->vx1,_X_->vy1,_X_->vx2,_X_->vy2,X,Y) / _X_->area;
			W1 = __AREA( _X_->vx2,_X_->vy2,_X_->vx0,_X_->vy0,X,Y) / _X_->area;
			W2 = __AREA( _X_->vx0,_X_->vy0,_X_->vx1,_X_->vy1,X,Y) / _X_->area;
			if(_X_->act){
				Z = (1/(W0 * (1/_X_->vz0) + W1 * (1/_X_->vz1) + W2 * (1/_X_->vz2)));
				U = (int)(( W0 * _X_->tx0 + W1 * _X_->tx1 + W2 * _X_->tx2 )*Z) % (_Model[_X_->Pmodel].FVM_HEAD.Num_Tam_T-1);
				V = (int)(( W0 * _X_->ty0 + W1 * _X_->ty1 + W2 * _X_->ty2 )*Z) % (_Model[_X_->Pmodel].FVM_HEAD.Num_Tam_T-1);
				*VIDEO = BilinealTex(U,V,_Model[_X_->Pmodel].FVM_HEAD.Num_Tam_T,_Model[_X_->Pmodel].TEX);
			}
			VIDEO++;
			DEPTH++;
			_X_->act = 0;
			_X_++;
		}
	}
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


int main(int argc, char const *argv[]){

	_GUB_NULL_FVM_(&_Model[0]);
	_GUB_NULL_FVM_(&_Model[1]);
	_GUB_NULL_FVM_(&_Model[2]);

	_GUB_LOAD_FVM_(&_Model[0],"extra.fvm");
	_GUB_LOAD_FVM_(&_Model[1],"soldier.fvm");
	_GUB_LOAD_FVM_(&_Model[2],"shambler.fvm");

	CreateCanvas(640,480);
	SetFillCanvas(0xffffffff);
	SetFillDepth(-999);

	float Camera0[3] = {1000,0,4000};
	float Camera1[3] = {-1000,0,4000};
	float Camera2[3] = {-300,300,1000};

	float Rot0[3] = {-170,0,0};
	float Rot1[3] = {-180,0,0};
	float Rot2[3] = {-180,0,0};

	DrawModel(0,Camera0,Rot0);
	DrawModel(1,Camera1,Rot1);
	DrawModel(2,Camera2,Rot2);

	SetShader();

	_GUB_FREE_FVM_(&_Model[0]);
	_GUB_FREE_FVM_(&_Model[1]);
	_GUB_FREE_FVM_(&_Model[2]);

	FlipBuff();
	FreeCanvas();

	return 0;
}