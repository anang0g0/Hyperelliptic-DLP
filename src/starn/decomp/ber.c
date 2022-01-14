#include "berlecamp.c"



//言わずもがな
int main(void)
{
    int i;
    unsigned short zz[N] = {0};
    OP f = {0}, r = {0}, w = {0}, r1 = {0};
    vec v, x = {0};
    MTX R = {0}, O = {0};
    unsigned short s[K + 1] = {0};

    if (K > N)
        printf("configuration error! K is bigger than N\n");

    // （謎）
    memset(mat, 0, sizeof(mat));

    // 公開鍵を生成する(Niederreiterとは異なる) // 鍵サイズK : Goppa Code
    R = pk_gen();
    // エラーベクトルの初期化
    memset(zz, 0, sizeof(zz));
    //重みTのエラーベクトルを生成する
    mkerr(zz, T);
    // 暗号文の生成(s=eH)
    x = sin3(zz, R);
    // 復号化１(m'=sS^{-1})
    r = dec(x.x);
    v = o2v(r);
    for (i = 0; i < K; i++)
        s[i + 1] = v.x[i];

    // Berlekamp-Massey Algorithm
    f = bma(s, K);
    x = chen(f);
    // 平文の表示(m=m'P^{-1})
    ero2(x);
    for(i=0;i<N;i++)
    if(zz[i]>0)
    printf("err=%d\n",i);
    //exit(1);
    wait();

    // debugging
    O = mk_pub(); // 鍵サイズ(K/2 Reed-Solomon)
    memset(zz, 0, sizeof(zz));
    //mkerr(zz, T);
    for(i=0;i<T;i++)
    zz[i]=1;
    r1 = sendrier2(zz, O);
    x = o2v(r1);
    for (i = 0; i < K; i++)
        s[i + 1] = x.x[i];
    //for (i = 0; i < K; i++)
    //    printf("%d,", s[i]);
    //printf("\n");
    f = bma(s, K);
    x = chen(f);
    ero2(x);
    printf("aaa\n");
    exit(1);
    for (i = 0; i < N; i++)
        if (zz[i] > 0)
            printf("%d,", i);
    printf("\n");

    if (odeg(f) < T)
    {
        printpol(o2v(r));
        printf("==goppa\n");
        for (i = 0; i < N; i++)
            printf("%d,", zz[i]);
        printf("\n");
        exit(1);
    }

    return 0;
}

/*
//言わずもがな
int main(void)
{
    int i;
    unsigned short zz[N] = {0};
    OP f = {0}, r = {0}, w = {0}, r1 = {0};
    vec v, x = {0};
    MTX R = {0}, O = {0};
    unsigned short s[K + 1] = {0};

    if (K > N)
        printf("configuration error! K is bigger than N\n");

    // （謎）
    memset(mat, 0, sizeof(mat));

    // 公開鍵を生成する(Niederreiterとは異なる) // 鍵サイズK : Goppa Code
    R = pk_gen();
    // エラーベクトルの初期化
    memset(zz, 0, sizeof(zz));
    //重みTのエラーベクトルを生成する
    mkerr(zz, T);
    // 暗号文の生成(s=eH)
    x = sin3(zz, R);
    // 復号化１(m'=sS^{-1})
    r = dec(x.x);
    v = o2v(r);
    for (i = 0; i < K; i++)
        s[i + 1] = v.x[i];

    // Berlekamp-Massey Algorithm
    f = bma(s, K);
    x = chen(f);
    // 平文の表示(m=m'P^{-1})
    ero2(x);
    for(i=0;i<N;i++)
    if(zz[i]>0)
    printf("err=%d\n",i);
    //exit(1);
    wait();

    // debugging
    O = mk_pub(); // 鍵サイズ(K/2 Reed-Solomon)
    memset(zz, 0, sizeof(zz));
    //mkerr(zz, T);
    for(i=0;i<T;i++)
    zz[i]=1;
    r1 = sendrier2(zz, O);
    x = o2v(r1);
    for (i = 0; i < K; i++)
        s[i + 1] = x.x[i];
    //for (i = 0; i < K; i++)
    //    printf("%d,", s[i]);
    //printf("\n");
    f = bma(s, K);
    x = chen(f);
    ero2(x);
    printf("aaa\n");
    exit(1);
    for (i = 0; i < N; i++)
        if (zz[i] > 0)
            printf("%d,", i);
    printf("\n");

    if (odeg(f) < T)
    {
        printpol(o2v(r));
        printf("==goppa\n");
        for (i = 0; i < N; i++)
            printf("%d,", zz[i]);
        printf("\n");
        exit(1);
    }

    return 0;
}
*/
