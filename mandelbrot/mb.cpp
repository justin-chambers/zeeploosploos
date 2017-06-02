#include <complex>
#include "PIMFuncs.h"
#include "Utilities.h"

void initialize_array(unsigned char ** arr);
int calculate_pixel_val(std::complex<float>& c);
void colorize(unsigned char ** raw,
              unsigned char ** red,
              unsigned char ** green,
              unsigned char ** blue,
              int m,
              int n);

int main(int argc, char ** argv) {
    // variables
    Utilities util;
    std::complex<float> pxl;
    int m, n;
    unsigned char ** img_r = NULL;
    unsigned char ** img_g = NULL;
    unsigned char ** img_b = NULL;

    // handle cmd ln
    util.handle_argc(argc, 4);
    m = util.get_positive_integer_from_cmd_line(argv[1]);
    n = util.get_positive_integer_from_cmd_line(argv[2]);

    // init array rows (n of them...)
    img_r = new unsigned char * [n];
    img_g = new unsigned char * [n];
    img_b = new unsigned char * [n];

    // init array row members (m of them...)
    for(int i = 0; i < n; ++i) {
        img_r[i] = new unsigned char[m];
        img_g[i] = new unsigned char[m];
        img_b[i] = new unsigned char[m];
    }
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            img_r[i][j] = 255;
            img_g[i][j] = 255;
            img_b[i][j] = 255;
        }
    }

    // calculate pixels
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            float re = (float) ((float(j) - float(m) / 2.0) * 4.0 / m);
            float im = (float) ((float(i) - float(n) / 2.0) * 4.0 / m);
            pxl.real(re);
            pxl.imag(im);
            img_r[i][j] = (const unsigned char)calculate_pixel_val(pxl);
        }
    }

    // colorize...
    colorize(img_r,img_r,img_g,img_b,m,n);

    // write img to file
    bool ret = pim_write_color(argv[3],
                               m,
                               n,
                               const_cast<const unsigned char **>(img_r),
                               const_cast<const unsigned char **>(img_g),
                               const_cast<const unsigned char **>(img_b));

     // clean up
    delete [] img_r;
    delete [] img_g;
    delete [] img_b;
    img_r = NULL;
    img_g = NULL;
    img_b = NULL;
    return 0;
}


int calculate_pixel_val(std::complex<float>& c){
    int pxl_val, max;
    std::complex<float> z;
    float tmp, length_squared;
    max = 256;
    z.real(0);
    z.imag(0);
    pxl_val = 0;
    do {
        tmp = (z.real()*z.real())-(z.imag()*z.imag())+c.real();
        z.imag(2*z.real()*z.imag()+c.imag());
        z.real(tmp);
        length_squared = z.real()*z.real()+z.imag()*z.imag();
        ++pxl_val;
    } while((length_squared < 4.0) && (pxl_val < max));
    if(pxl_val < max) return pxl_val;
    else return 0;
}

void colorize(unsigned char ** raw,
              unsigned char ** red,
              unsigned char ** green,
              unsigned char ** blue,
              int m,
              int n) {
    int r = 54, g = 73, b = 104;
    for(int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            int value = (int)raw[i][j];
            int tmp_r=0, tmp_g=0, tmp_b=0;
            if(value != 0) {
                tmp_r = abs(r-value);
                tmp_g = abs(g-value);
                tmp_b = abs(b-value);
            }
            red[i][j] = (const unsigned char)(tmp_r);
            green[i][j] = (const unsigned char)(tmp_g);
            blue[i][j] = (const unsigned char)(tmp_b);
        }
    }
}

//everything i know i learned in kindergarten
//who moved my cheese