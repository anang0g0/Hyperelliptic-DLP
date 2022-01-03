# Hyperelliptic-DLP (･∀･) ﾔｺﾋﾞﾔｰﾝ!!

# 20220103

ZZ P = to_ZZ("2923003274661805836407369665432566039311865180529"); //基礎となる有限体  
ZZ aa = to_ZZ("371293"); // y^2=x^5+aax  
ZZ J =to_ZZ("8543948143683640329580084318401338115672828124663448275867130387651937373152534160174163969676194");  //ヤコビアンの位数  

このように設定してあります。  
ヤコビアンの位数は、

3^2 * 5 * 2393 * 770767 * 1387123 * 25620149 * 2896567510666771668598410666311314487756987819570261851119062333047490183

このように因数分解できるので、この中の最大素数は241bitあります。（間違っているかもしれません）  
つまり、この暗号は240ビット以上の安全性を持つことがわかります。

ここまで短期間で実装できたのは、私は膜リースのために2年間費やして作った、多項式ライブラリがあったおかげで、このように多項式ライブラリを応用に使えるという証拠になったわです。
そして今年の目標は多変数多項式ライブラリにしようと思います。

ちなみにここで実装しているのはカントールアルゴリズムです。

コンパイルにはNTL6.2を使っています。
先のNTLをビルドしてインストールして使ってください。


https://eprint.iacr.org/2002/181.pdf

https://www.researchgate.net/publication/267120817_Fast_Cryptography_in_Genus_2

なんか256ビットくらいの超楕円曲線でも計算できてるみたいなので一応完成ということに。
次何やろうかな。

# 20220102

家の中に引きこもっているので正月という気分になりません。
一日中デバッグしてます。

なんだかバグを取るのが趣味のようになってしまい、こうなるとプログラミングというものも中毒性があります。

前は動かないとイライラしたり、やたらと不安になったりしましたが、今は余裕を持ってバグ取りができます。
ある程度上達したからでしょうか？

今はランダム因子の生成のバグをとっています。
因子の多項式の次数を減らすための処理で、逆に多項式の次数が上がってしまうなど、数え切れバイバグがあります。

正月もパソコンで終わってしまうのだろうか。

# 20220101

あけましておめでとうございます。
あまり正月って感じではないですが、今年も健康に過ごせたらいいなと思っています。

巨大整数バージョンに対応しました。
目下のところランダムに因子をとってこようとすると、確率アルゴリズムでは膨大な時間がかかってしまうので、それを改善しています。

ドキュメントも書きたいですが、巨大整数までできたので今日はおしまい。

美味しいおしるこに、お雑煮。

食い物で正月を演出しています。

# 20211231

次は点数えなんかもやってみたいです。

https://eprint.iacr.org/2011/306.pdf

なんか出来ちゃった感があるのでお正月はゆっくりしたいです。（仕事でもないのに）

cantorアルゴリズムの加法と2倍点計算が（ゆるく）できるようになりました。
まだ加法計算にバグがあるのでそれは直さないといけないです。
でも元のプログラムがかなりアバウトだったせいで、細かいバグを修正して、かなり精度は上がったはずです。

素数を法とする有限体上の多項式の計算なので、超楕円だけでなくWild McElieceなんかにも使えたりして意外と応用範囲が広めで便利です。

2倍点の計算はチェコの論文を参考にしました。

https://www.vut.cz/www_base/zav_prace_soubor_verejne.php?file_id=30949

ドキュメントはこれから書きます。

# 20211230

まだデバッグの段階なのであまり細かくは書きませんが、プログラム本体はcantorディレクトリのjacobian.cです。

# 20211224

バイナリのときにはうまく行っていた逆多項式の計算が素体に変えたら動かなくなったので、逆多項式を計算するなにかいい方法がないかと調べたら拡張ユークリッドが使えるというので一日かけてバグをとった。

そのおかげでカントールアルゴリズムを使って超楕円曲線のヤコビアンの加法を実装できるようになった。（数学的な理解ではなく計算だけ）
超楕円の場合検証用のデータがないのであっているかどうかがわかりにくい。

結局最後は自己解決するしないのか。同じ場所に座り続けているので畳が擦れてきたので座布団を敷いて作業しているｗ

# 20211223

楕円ができたんだから超楕円なんか楽勝！って言いたいところだが、素体上に自分の書いたプログラムを変更しようとすると途端に動かなくなる。
今はもう前のようなエネルギーが枯渇しているので、あんなに頑張れない。あれはやはり躁状態だったからできたんだと思う。
仕方がないので誰かがPythonで書いたやつを使おうとしてみたり、sagemathのスクリプトで我慢しようとしてみたりもしたが、いずれも成功せず。一体どうすればいいのか。

年末から来年にかけては、自分の多項式演算プログラムに本格的な手直しを入れるかもしれない。（気合があれば）

# 20211221

なんでも量子計算機の脅威は当分ないとかで、来年以降の暗号の主軸は多機能暗号と軽量暗号になるらしい。

その他別枠でPQCとのこと。

従来方式ならぜひ超楕円をやりたい。
まだデバックしてるけどな。


# 20211220

このまま超楕円曲線の位数計算まで行けば面白いかも。

暗号としてはすでに価値がないとしても、昔わからなかったﾔｺﾋﾞﾔｰﾝの無念を晴らすことができそう。

# 

種数２と３しかできてません。

https://scholar.uwindsor.ca/cgi/viewcontent.cgi?article=6719&context=etd

４までできますがやりません。
