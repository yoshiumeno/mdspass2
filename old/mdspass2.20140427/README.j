mdspass2 インストールマニュアル

１）freeglut のインストール
freeglut と freeglut-devel のパッケージがインストールされている必要があります。
rootになって
# yum install freeglut
# yum install freeglut-devel
でＯＫのはず。

２）glui のインストール
glui は http://glui.sourceforge.net/ より入手可能。
glui-2.36.tgz をダウンロード、コンパイル。
　srcディレクトリに行き、makefileを修正
　(2)の FreeGLUTのところのコメントアウトを外し、
    /usr/X11R6/lib -lfreeglut を /usr/lib -lglut に、
    /usr/X11R6/include を /usr/include に直す。
  (3)のところはコメントアウト。
make し、
所定の場所に libglui.a と glui.h をコピー。
私の環境では、
 /usr/lib/libglui.a
 /usr/include/GL/glui.h

＝＝＝Linux版GLUIのバグとその修正について＝＝＝
（面倒な方は「バグその１・その２」は読まずに「バグその３」へジャンプしてください）
◇GLUIのバグその１◇
システムによっては，ファイルブラウザでダブルクリックが効かないというバグに
遭遇する（example6で，左下のボックスに現れるファイルをダブルクリックしても
右側の窓にファイル内容が表示されない）．
その際は，
src/glui_list.cpp の126行め
  if (last_line == curr_line && (ms - last_click_time) < 300) {
を
  if (last_line == curr_line && (int)(ms - last_click_time) < 300) {
のように修正する必要がある．

◇GLUIのバグその２◇
・ファイルブラウザでディレクトリ変更ができない，
・一番トップに表示されるファイルが選択できない，
というバグに遭遇する（上記と同様，example6でテストされたい）．
その際は，
src/glui_filebrowser.cppの79行め
  if (this_item > 0) { /* file or directory selected */
を
  if (this_item >= 0) { /* file or directory selected */
に，81行め
    if (selected[0] == '/' || selected[0] == '\\') {
を
    if (selected[strlen(selected)-1] == '/' || selected[0] == '\\') {
に，84行め
        chdir(selected+1);
を
        chdir(selected);
に修正する必要がある．

◇GLUIのバグその３◇
ファイルブラウザでファイル名がアルファベット順に表示されず，出鱈目な順序なので
非常に使いにくい．これを修正するには新たな関数を組込むなど比較的大きな改修が
必要であったので，ここでは書ききれない．
mdspass2のディレクトリ以下に
  glui-2.36_bugfix2/
というディレクトリがあるが，そこに修正版のglui一式を入れておいたので，
そこに移動してから上記２）の手順を実行すればよい．

３）libpngのインストール
libpng と libpng-devel がインストールされている必要があります．
yum install libpng-devel で完了．

注）makeの際、環境によっては、-lXmu, -lXiなどがないと言われる。
 その場合には例えば
 yum install libXmu-devel, yum install libXi-devel などでインストール。

４）lapackのインストール
rootになって
# yum install lapack
# yum install blas
# yum install lapack-devel
とする．

５）CUDA版のビルド（開発中なので読み飛ばしてください）
Makefileの最初
CC = g++　をコメントアウトし，　CC = nvcc -DCUDA　のコメントアウトを外す．
その後，
make clean; make cuda
を実行する．
