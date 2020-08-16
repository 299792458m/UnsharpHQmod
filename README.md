# UnsharpHQmod
a fork of list's UnsharpHQ  
  UnsharpHQ is a video plugin for Avisynth

## PARAMETERS
(All parameters are optional).
* THRESHOLD: (Default = 20)  The value to determine whether or not to sharp a pixel based on the luminance change between their neighbors. Any value >0 can be used.
* SHARPSTR : (Default = 4.0) The sharp strength for the pixels to be sharped. Values from 2 to 10 are recommended.
* SMOOTH   : (Default = 0.5) The smooth strength to pixels not sharped. Use this with 0 if you don't want any smooth, or just leave it with default.

	EXAMPLES
Unsharphq()           #All by default, this is useful in most cases.
Unsharphq(20,15.0)    #Very strong sharp.
--------------------------------------------------------------------------
### modのメモ
* MODE	:0 default(従来相当)  
		:1 mod版で追加 パラメータの効果が変わるので分けてる
* SHOW	:いまいち使い方がわからないのでHBDではちゃんとした対応をしていない opt=1の時のみ対応
* OPT	:最適化 0=auto 1=C 2:SSE2 3:AVX2
		
### 変更履歴 mod by 299792458m
- mod 200816
  - Avisynth2.6+化
  - HBD(8-16 YUVのみ)対応＋AVX2最適化、パラメータにopt追加
  - asmのintrinsic化
- mod 190119
  - ループの範囲を修正
  - 補正重みを4→8に(これによりmode0でも出力が微妙に変更される)
  - Mode1での補正重みを中央点と8近傍の差分がTHRESHOLDの時に1になるように修正
- mod 181224	いろいろ

