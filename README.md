# Hyperelliptic-DLP

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
