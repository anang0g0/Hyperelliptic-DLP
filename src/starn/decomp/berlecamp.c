#include "common.c"

OP mkd(OP w, int kk)
{
    int i, j, k, l, ii = 0;

    unsigned short tr[N] = {0};
    unsigned short ta[N] = {0};
    vec v = {0};
    unsigned short po[K + 1] = {1, 0, 1, 0, 5};
    //OP w={0};
    OP r = {0};

aa:

    //printf("\n");
    memset(mat, 0, sizeof(mat));
    //既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
    //既約多項式しか使わない。

    l = -1;
    ii = 0;
    // irreducible goppa code (既役多項式が必要なら、ここのコメントを外すこと。)
    
    while (l == -1)
    {
        w = mkpol();
        l = ben_or(w);
        printf("irr=%d\n", l);
        if (ii > 300)
        {
            printf("too many error\n");
            exit(1);
        }
        ii++;
        //
    }

    // separable goppa code
    //w = mkpol();
    r = w;
    //  r=omul(w,w);
    memset(ta, 0, sizeof(ta));
    //w = setpol(g, K + 1);
    printpol(o2v(r));
    //printf(" =poly\n");

    //多項式の値が0でないことを確認
    for (i = 0; i < N; i++)
    {
        ta[i] = trace(r, i);
        if (ta[i] == 0)
        {
            printf("trace 0 @ %d\n", i);
            //fail = 1;
            goto aa;
        }
    }
    for (i = 0; i < N; i++)
    {
        tr[i] = oinv(ta[i]);
        //printf("%d,", tr[i]);
    }
    memset(g, 0, sizeof(g));
    g[0] = 1;

    //多項式を固定したい場合コメントアウトする。
    oprintpol(r);
    printf("\n");
    printsage(o2v(r));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");
    memset(v.x, 0, sizeof(v.x));
    //  v=rev(w);
    van(kk);
    //  v=o2v(w);
    ogt(g, kk);

    //wait();

    //#pragma omp parallel for

    printf("\nすげ、オレもうイキそ・・・\n");
    //keygen(g);
    //exit(1);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < kk; j++)
        {
            mat[i][j] = vb[j][i];
        }
    }

    //printf("\n");
    //exit(1);
    /*
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < kk; i++)
            printf("%d,", mat[j][i]);
        printf("\n");
    }
    //exit(1);
    //wait();
*/

    return w;
}


unsigned short HH[N][K];
unsigned short TE[N][K / 2 + 1];

MTX toByte2(MTX SH, int kk)
{
    vec v = {0};
    int i, j, k, cnt;
    MTX R = {0};

    memset(HH, 0, sizeof(HH));
    printf("HH=");
    //exit(1);
    for (i = 0; i < N; i++)
    {
        printf("%d\n",i);
        //#pragma omp parallel for
        for (j = 0; j < kk; j++)
        {
            cnt = 0;
            for (k = j * E; k < j * E + E; k++)
                v.x[cnt++] = SH.x[i][k];

            HH[i][j] = v2i(v);
            R.x[i][j] = v2i(v);
            //printf("%d,", HH[i][j]);
            //= BH[j][i];
        }
        //fwrite(dd, 1, E * K, ff);
        //printf("\n");
    }
    printf("end of byte\n");
    //exit(1);
    //wait();

    return R;
}

MTX bd2()
{
    int i, j, k, l;
    unsigned char dd[E * K] = {0};
    FILE *ff;
    vec v = {0};
    MTX R = {0};

    //ff = fopen("Hb.key", "wb");

    //memset(BB.z,0,sizeof(BB.z));
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K / 2 + 1; j++)
        {
            l = bm[i][j];
            printf("bm==%d %d\n", l, j);

            v = i2v(l);
            //#pragma omp parallel for
            for (k = 0; k < E; k++)
            {
                R.x[i][j * E + k] = v.x[k];
                //l = (l >> 1);
            }
        }
    }

    return R;
}


//Niederreiter暗号の公開鍵を作る(RS)
MTX mk_pub()
{
    int i, j, k, l;
    FILE *fp;
    unsigned char dd[E * K] = {0};
    OP w = {0};
    MTX Z = {0}, FX = {0}, G_bin = {0}, G_int = {0}, R = {0}, O = {0};

    w = mkd(w, K * 2);
    //w = mkg(K);
    half(K / 2 + 1);

    oprintpol(w);
    printf("\n");
    printsage(o2v(w));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");

    R = bd2();
    printf("deco_rev= ");
    for (i = 0; i < (K / 2 + 1) * E; i++)
        printf("%d,", R.x[1][i] ^ R.x[2][i] ^ R.x[3][i]);
    printf("\n");
    Pgen();
    do
    {
        memset(inv_S.x, 0, sizeof(inv_S.x));
        memset(S.x, 0, sizeof(S.x));
        for (i = 0; i < (K/2+1) * E; i++)
        {
            for (j = 0; j < (K/2+1) * E; j++)
                S.x[i][j] = xor128() % 2;
        }
    } while (mkRS(S, inv_S.x) == -1);
    //mkS();
    //O = toByte(R, K / 2 + 1);
    //  exit(1);
    Z = mulmat(S, R, 2);
    printf("Zz=\n");
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < (K / 2 + 1) * E; i++)
        {
            //G.z[j][i] = Z.w[P[j]][i];
            G_bin.x[j][i] = Z.x[P[j]][i];
            //Z.x[j][i]=Z.x[j][i];
            //printf("%d.",inv_S.w[j][i]);
        }
        //printf("\n");
    }
    //printf("\n");


    FX = toByte2(G_bin, K / 2 + 1);
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K / 2 + 1; j++)
            G_int.x[i][j] = FX.x[i][j];
    }

    return G_int;
}


//バイナリ型パリティチェック行列を生成する
MTX bdet2()
{
    int i, j, k, l;
    unsigned char dd[E * K] = {0};
    FILE *ff;
    MTX R = {0};

    //ff = fopen("Hb.key", "wb");

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
        {
            l = mat[i][j];
            //#pragma omp parallel for
            for (k = 0; k < E; k++)
            {
                R.x[i][j * E + k] = l % 2;
                l = (l >> 1);
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        //#pragma omp parallel for
        for (j = 0; j < E * K; j++)
        {
            printf("%d,", R.x[i][j]);
            //dd[j] = BH[j][i];
        }
        //fwrite(dd, 1, E * K, ff);
        printf("\n");
    }

    //fclose(ff);
    return R;
}


//Niederreiter暗号の公開鍵を作る(Goppa)
MTX pk_gen()
{
    int i, j, k, l;
    FILE *fp;
    unsigned char dd[E * K] = {0};
    OP w = {0};
    MTX R = {0}, R_bin = {0}, O = {0}, Q = {0}, O_bin = {0};

    w = mkd(w, K);
    //w = mkg(K);
    //half(K / 2 + 1);

    oprintpol(w);
    printf("\n");
    printsage(o2v(w));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");

    Q = bdet2();

    Pgen();
    do
    {
        memset(inv_S.x, 0, sizeof(inv_S.x));
        memset(S.x, 0, sizeof(S.x));
        for (i = 0; i < K * E; i++)
        {
            for (j = 0; j < K * E; j++)
                S.x[i][j] = xor128() % 2;
        }
    } while (is_reg(S, inv_S.x) == -1);

    //makeS();
    //  exit(1);
    H = mulmat(S, Q, 1);
    for (i = 0; i < K * E; i++)
    {
        for (j = 0; j < N; j++)
        {
            //O.z[j][i] = H.z[P[j]][i];
            O_bin.x[j][i] = H.x[P[j]][i];
        }
    }
    R = toByte2(O_bin, K);
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K * E; j++)
        {
            R_bin.x[i][j] = O_bin.x[i][j];
            //printf("%d,",R_bin.x[i][j]);
        }
    }

    //exit(1);

    return R;
}

void half(int kk)
{
    int i, j, k;

    for (i = 0; i < N; i++)
        bm[i][0] = mat[i][0];
    for (i = 1; i < kk; i++)
    {
        for (j = 0; j < N; j++)
            bm[j][i] = mat[j][i * 2 - 1];
    }
    /*
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < kk+1; j++)
            printf("%d,", bm[i][j]);
        printf("  ==bm\n");
    }
*/
    //exit(1);
}

OP bma(unsigned short s[], int kk)
{
    int i, j, k, ll = 0, l, d[2 * K + 1] = {0};
    OP lo[2 * K + 1] = {0}, b[2 * K + 1] = {0}, t[2 * K + 1] = {0}, a = {0}, f = {0}, h = {0}, g = {0}, hh = {0};
    vec v = {0}, x = {0}, w = {0};

    x.x[1] = 1;
    h = v2o(x);
    v.x[0] = 1;
    f = v2o(x);
    lo[0] = v2o(v);
    b[0] = lo[0];
    ll = 0;
    for (j = 1; j < T * 2 + 1; j++)
    {
        v = o2v(lo[j - 1]);
        k = 0;
        //printpol(v);
        //printf(" ==lo\n");

        l = deg(o2v(lo[j - 1]));
        for (i = 1; i < l + 1; i++)
        {
            k ^= gf[mlt(fg[v.x[i]], fg[s[j - i]])];
            //printf("v[%d]=%d\n", i, v.x[i]);
        }
        d[j] = s[j] ^ k;
        //printf("d[%d]=%d\n", j, d[j]);
        if (d[j] == 0)
        {
            lo[j] = lo[j - 1];
            b[j] = omul(b[j - 1], h);
            //ll=j-1;
        }
        else //if (d[j] > 0)
        {
            g = omul(kof(d[j], h), b[j - 1]);
            t[j] = oadd(lo[j - 1], g);
            lo[j] = t[j];
            if (ll * 2 > (j - 1))
            {
                //lo[j]=t[j];
                b[j] = omul(b[j - 1], h);
            }
            else //if(2*ll <= j)
            {
                //printpol(o2v(t[j]));
                //printf("==t[%d]\n", j);
                b[j] = kof(gf[oinv(d[j])], lo[j - 1]);
                //lo[j]=t[j];
                ll = j - ll;

                if (j == 2 * T)
                {
                    if (!(d[T * 2 - 1] == 0 && d[T * 2 - 3] == 0 && odeg(lo[j - 1]) == T) || !(odeg(lo[j - 1]) == T))
                    {
                        if ((d[T * 2 - 1] == 0 && odeg(lo[j - 2]) == T - 1))
                        {
                            lo[j - 1] = omul(lo[j - 2], h);
                            //printpol(o2v(lo[j - 1]));
                            //printf("\n");
                        }
                    }
                    break;
                }
            }
        }
        printf("l=%d\n", ll);
        k = 0;
        //printpol(o2v(b[j]));
        //printf(" ==b[%d]\n", j);
    }

    k = 0;
    int count = 0;
    //printpol(o2v(lo[j - 1]));
    //printf(" ==coef\n");
    if (odeg(lo[j - 1]) == T)
    {
        x = chen(lo[j - 1]);
    }
    else
    {
        printf("baka***\n");
        exit(1);
        //return -1;
    }
    //exit(1);
    for (i = 0; i < deg(x) + 1; i++)
    {
        if (x.x[i] >= 0)
        {
            printf("xx[%d]=1\n", (x.x[i]));
            count++;
        }
        //

        if (x.x[i] == 0)
            k++;
        if (k > 1)
        {
            printf("baka0\n");
            //printvec((x));
            //for (i = 0; i < N; i++)
            //printf("%d,", zz[i]);
            exit(1);
            //return f;
        }
    }
    if (count < T)
    {
        printf("vaka in bms %d\n", count);
        exit(1);
    }

    //return count;
    return lo[j - 1];
}


vec newhalf(unsigned short e[K])
{
    int i, j, k;
    vec v = {0};
    unsigned short t[K]={0};

    for (i = 0; i < K / 2 + 1; i++){
    t[i]=e[i];
        printf("e=%d\n", e[i]);
    }
    for (i = 0; i < K / 2 + 1; i++){
        printf("t=%d\n", t[i]);
    }
    //exit(1);

    v.x[0] = t[0];
    v.x[1] = t[1];
    k = 2;
    for (i = 2; i < K; i++)
    {
        if (i % 2 == 1)
        {
            v.x[i] = t[k];
            k++;
        }
        if (i % 2 == 0)
            v.x[i] = gf[mlt(fg[v.x[i / 2]], fg[v.x[i / 2]])];
    }

    return v;
}

vec bfd(unsigned short ss[])
{
    int i, j, k, count = 0;
    vec v = {0};
    OP s = {0};
    unsigned ch[K * E * 2] = {0};
    unsigned char h2o[K * E * 2] = {0};

    //count=(K/2+1)*E-1;
    for (i = 0; i < (K / 2) + 1; i++)
    {
        v = i2v(ss[i]);
        for (j = 0; j < E; j++)
        {
            ch[count] = v.x[j];
            count++;
        }
    }
    printf("bfd_bin=\n");
    for (i = 0; i < (K / 2 + 1) * E; i++)
        printf("%d", ch[i]);
    printf("\n");
    //exit(1);

    unsigned short uk[K] = {0};

    for (i = 0; i < (K / 2 + 1) * E; i++)
    {
        for (j = 0; j < (K / 2 + 1) * E; j++)
            h2o[i] ^= (ch[j] & inv_S.x[i][j]);
    }

    printf("deco_bin=\n");
    for (i = 0; i < (K / 2 + 1) * E; i++)
        printf("%d,", h2o[i]);
    printf("\n");

    //count=(K/2+1)*E-1;

    count = 0;
    for (i = 0; i < (K / 2) + 1; i++)
    {
        memset(v.x, 0, sizeof(v.x));
        for (j = 0; j < E; j++)
        {
            v.x[j] = h2o[count];
            count++;
        }
        uk[i] = v2i(v);
    }

    printf("bm_int=\n");
    for (i = 0; i < K / 2 + 1; i++)
        printf("%d,", bm[1][i] ^ bm[2][i] ^ bm[3][i]);
    printf("\n");

    printf("u_int=\n");
    for (i = 0; i < K / 2 + 1; i++)
        printf("%d,", uk[i]);
    printf("\n");
    //exit(1);

    s = setpol(uk, K / 2 + 1);
    v = o2v(s);

    return v;
}



OP sendrier2(unsigned short zz[N], MTX L)
{
    unsigned short syn[K + 1] = {0}, s[K + 1] = {0}, rt[K / 2 + 1] = {0}, uu[(K / 2 + 1) * E] = {0}, es[(K / 2 + 1) * E] = {0};
    int i, j, k, count = 0;
    OP f = {0}, w = {0};
    vec v = {0}, x = {0}, u = {0}, t = {0};
    unsigned short tmp[(K / 2 + 1) * E] = {0}, m[K + 1] = {0};

    for (j = 0; j < N; j++)
    {
        if (zz[j] > 0)
        {
            //for(i=0;i<(K/2)+1;i++){

            //memcpy(syn, L.w[j], sizeof(syn));

            for (k = 0; k < K / 2 + 1; k++)
            {
                rt[k] ^= L.x[j][k];
                //rt[k] = bm[j][k];
            }
        }
    }
    for(i=0;i<K/2+1;i++)
    printf("%d,",rt[i]);
    printf("\n");
    //exit(1);

    x = bfd(rt);
    for (i = 0; i < K / 2 + 1; i++){
        u.x[K / 2 - i] = x.x[i];

    }
    for(i=0;i<K/2+1;i++)
        printf("%d,%d\n",u.x[i],x.x[i]);
    printf("\n");
    //exit(1);

    printf("rt=\n");
    for (i = 0; i < K / 2 + 1; i++)
        printf("%d,", u.x[i]);
    printf("\n");
    printf("bm_in se2 == %d || ", j);
    //exit(1);
    v = newhalf(u.x);

    printf("P= ");
    for (i = 0; i < N; i++)
        printf("%d,", P[i]);
    printf("\n");

    memset(s,0,sizeof(s));
    for (i = 0; i < K; i++)
        s[i + 1] = v.x[i];
    
    printf("rt_deco= ");
    f=bma(s, K);
    x=chen(f);
    ero2(x);
    //exit(1);
    wait();


    f = setpol(v.x, K);

    return f;
}

vec sin3(unsigned short zz[], MTX R)
{
    int i, j;
    OP s = {0};
    vec v = {0};
    unsigned short ss[K] = {0};

    for (i = 0; i < N; i++)
    {
        if (zz[i] > 0)
        {
            for (j = 0; j < K; j++)
            {
                //ss[j]
                v.x[j] ^= R.x[i][j];
                printf("%d,", R.x[i][j]);
            }
        }
        printf("\n");
    }

    return v;
}

