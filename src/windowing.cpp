#include <plugin.h>

typedef struct _fft {
  OPDS h;
  ARRAYDAT *out;
  ARRAYDAT *in, *in2;
  MYFLT *f;
  MYFLT b;
  int32_t n;
  void *setup;
  AUXCH mem;
} FFT;


static inline void tabinit(CSOUND *csound, ARRAYDAT *p, int size)
{
    size_t ss;
    if (p->dimensions==0) {
        p->dimensions = 1;
        p->sizes = (int32_t*)csound->Malloc(csound, sizeof(int32_t));
    }
    if (p->data == NULL) {
        CS_VARIABLE* var = p->arrayType->createVariable(csound, NULL);
        p->arrayMemberSize = var->memBlockSize;
        ss = p->arrayMemberSize*size;
        p->data = (MYFLT*)csound->Calloc(csound, ss);
        p->allocated = ss;
    } else if( (ss = p->arrayMemberSize*size) > p->allocated) {
        p->data = (MYFLT*) csound->ReAlloc(csound, p->data, ss);
        p->allocated = ss;
    }
    if (p->dimensions==1) p->sizes[0] = size;
    //p->dimensions = 1;
}

int32_t shiftin_init(CSOUND *csound, FFT *p) {
    int32_t sizs = CS_KSMPS;
    if(p->out->sizes[0] < sizs)
       tabinit(csound, p->out, sizs);
    p->n = 0;
    return OK;
}

int32_t shiftin_perf(CSOUND *csound, FFT *p) {
   IGN(csound);
    uint32_t  siz =  p->out->sizes[0], n = p->n;
    MYFLT *in = ((MYFLT *) p->in);
    if (n + CS_KSMPS < siz) {
      memcpy(p->out->data+n,in,CS_KSMPS*sizeof(MYFLT));
    }
    else {
      int32_t num = siz - n;
      memcpy(p->out->data+n,in,num*sizeof(MYFLT));
      memcpy(p->out->data,in+num,(CS_KSMPS-num)*sizeof(MYFLT));
    }
    p->n = (n + CS_KSMPS)%siz;
    return OK;
}


int32_t shiftout_init(CSOUND *csound, FFT *p) {
    int32_t siz = p->in->sizes[0];
    p->n = ((int32_t)*((MYFLT *)p->in2) % siz);
    if (UNLIKELY((uint32_t) siz < CS_KSMPS))
      return csound->InitError(csound, "%s", Str("input array too small\n"));
    return OK;
}

int32_t shiftout_perf(CSOUND *csound, FFT *p) {
    IGN(csound);
    uint32_t siz =  p->in->sizes[0], n = p->n;
    MYFLT *out = ((MYFLT *) p->out);

    if (n + CS_KSMPS < siz) {
      memcpy(out,p->in->data+n,CS_KSMPS*sizeof(MYFLT));
    }
    else {
      int32_t num = siz - n;
      memcpy(out,p->in->data+n,num*sizeof(MYFLT));
      memcpy(out+num,p->in->data,(CS_KSMPS-num)*sizeof(MYFLT));
    }
    p->n = (n + CS_KSMPS)%siz;
    return OK;
}


struct windowing : csnd::Plugin<1, 1> {
    static constexpr char const *otypes = "k[]";
    static constexpr char const *itypes = "k[]";
    FWT<MYFLT> fwt;
    
    int init() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        out.init(csound, in.len());
        return OK;
    }
    
    int kperf() {
        csnd::Vector<MYFLT> &out = outargs.vector_data<MYFLT>(0);
        csnd::Vector<MYFLT> &in = inargs.vector_data<MYFLT>(0);
        std::copy(in.begin(), in.end(), out.begin());
        fwt.fwtanal(out);
        return OK;
    }
};

int32_t init_window(CSOUND *csound, FFT *p) {
    int32_t   N = p->in->sizes[0];
    int32_t   i,type = (int32_t) *p->f;
    MYFLT *w;
    tabinit(csound, p->out, N);
    if (p->mem.auxp == 0 || p->mem.size < N*sizeof(MYFLT))
      csound->AuxAlloc(csound, N*sizeof(MYFLT), &p->mem);
    w = (MYFLT *) p->mem.auxp;
    switch(type) {
    case 0:
      for (i=0; i<N; i++) w[i] = 0.54 - 0.46*cos(i*TWOPI/N);
      break;
    case 1:
    default:
      for (i=0; i<N; i++) w[i] = 0.5 - 0.5*cos(i*TWOPI/N);
    }
    return OK;
}

int32_t perf_window(CSOUND *csound, FFT *p) {
    IGN(csound);
    int32_t i,end = p->out->sizes[0], off = *((MYFLT *)p->in2);
    MYFLT *in, *out, *w;
    in = p->in->data;
    out = p->out->data;
    w = (MYFLT *) p->mem.auxp;
    /*while (off < 0) off += end;
     for (i=0;i<end;i++)
      out[(i+off)%end] = in[i]*w[i];*/
    if(off) off = end - off;
    for(i=0;i<end;i++)
      out[i] = in[i]*w[(i+off)%end];
    return OK;

}



static inline void tabensure2D(CSOUND *csound, ARRAYDAT *p,
                               int32_t rows, int32_t columns)
{
    if (p->data==NULL || p->dimensions == 0 ||
        (p->dimensions==2 && (p->sizes[0] < rows || p->sizes[1] < columns))) {
      size_t ss;
      if (p->data == NULL) {
        CS_VARIABLE* var = p->arrayType->createVariable(csound, NULL);
        p->arrayMemberSize = var->memBlockSize;
      }
      ss = p->arrayMemberSize*rows*columns;
      if (p->data==NULL) {
        p->data = (MYFLT*)csound->Calloc(csound, ss);
        p->dimensions = 2;
        p->sizes = (int32_t*)csound->Malloc(csound, sizeof(int32_t)*2);
      }
      else p->data = (MYFLT*) csound->ReAlloc(csound, p->data, ss);
      p->sizes[0] = rows;  p->sizes[1] = columns;
    }
}
// setrow
int32_t set_rows_init(CSOUND *csound, FFT *p) {
    int32_t sizs = p->in->sizes[0];
    int32_t row = *((MYFLT *)p->in2);
    tabensure2D(csound, p->out, row+1, sizs);
    return OK;
}

int32_t set_rows_perf(CSOUND *csound, FFT *p) {
    int32_t start = *((MYFLT *)p->in2);
    if (UNLIKELY(start < 0 || start >= p->out->sizes[0]))
      return csound->PerfError(csound, &(p->h),
                                 Str("Error: index out of range\n"));
    int32_t bytes =  p->in->sizes[0]*sizeof(MYFLT);
    start *= p->out->sizes[1];
    memcpy(p->out->data+start,p->in->data,bytes);
    return OK;
}


// getrow
int32_t rows_init(CSOUND *csound, FFT *p) {
    if (p->in->dimensions == 2) {
      int32_t siz = p->in->sizes[1];
      tabinit(csound, p->out, siz);
      return OK;
    }
    else
      return csound->InitError(csound, "%s",
                               Str("in array not 2-dimensional\n"));
}

int32_t rows_perf(CSOUND *csound, FFT *p) {
    int32_t start = *((MYFLT *)p->in2);
    if (LIKELY(start < p->in->sizes[0])) {
      int32_t bytes =  p->in->sizes[1]*sizeof(MYFLT);
      start *= p->in->sizes[1];
      memcpy(p->out->data,p->in->data+start,bytes);
      return OK;
    }
    else return csound->PerfError(csound,  &(p->h),
                                  Str("requested row is out of range\n"));
}