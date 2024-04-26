## direct copy paste of ipa table from here:  https://en.wikipedia.org/wiki/Burmese_alphabet
library(tidyverse)
## gc@(aspirated, unaspirated, voiced, nasal) -> frame@{gc$aspirated, ...}
## Unaspirated (သိထိလ)	Aspirated (ဓနိတ)	Voiced (လဟု)	Nasal (နိဂ္ဂဟိတ)
tribble(~group_name, ~bgn, ~aspirated, ~bname, ~voiced, ~nasal,
"Velars","ကဏ္ဍဇ",
ကဝဂ်	က	k	/k/	ခ	hk	/kʰ/	ဂ	g	/ɡ/	ဃ	gh	/ɡ/	င	ng	/ŋ/
ကကြီး [ka̰ dʑí]	ခကွေး [kʰa̰ ɡwé]	ဂငယ် [ɡa̰ ŋɛ̀]	ဃကြီး [ɡa̰ dʑí]	င [ŋa̰]
Palatals
(တာလုဇ)
စဝဂ်	စ	c	/s/	ဆ	hc	/sʰ/	ဇ	j	/z/	ဈ	jh	/z/	ဉ / ည	ny	/ɲ/
စလုံး [sa̰ lóʊɰ̃]	ဆလိမ် [sʰa̰ lèɪɰ̃]	ဇကွဲ [za̰ ɡwɛ́]	ဈမျဉ်းဆွဲ [za̰ mjɪ̀ɰ̃ zwɛ́]	ညကလေး/ ညကြီး [ɲa̰ dʑí]
Alveolars
(မုဒ္ဓဇ)
ဋဝဂ်	ဋ	t	/t/	ဌ	ht	/tʰ/	ဍ	d	/d/	ဎ	dh	/d/	ဏ	n	/n/
ဋသန်လျင်းချိတ် [ta̰ təlɪ́ɰ̃ dʑeɪʔ]	ဌဝမ်းဘဲ [tʰa̰ wʊ́ɰ̃ bɛ́]	ဍရင်ကောက် [da̰ jɪ̀ɰ̃ ɡaʊʔ]	ဎရေ မှုတ် [da̰ jè m̥oʊʔ]	ဏကြီး [na̰ dʑí]
Dentals
(ဒန္တဇ)
တဝဂ်	တ	t	/t/	ထ	ht	/tʰ/	ဒ	d	/d/	ဓ	dh	/d/	န	n	/n/
တဝမ်းပူ [ta̰ wʊ́ɰ̃ bù]	ထဆင်ထူး [tʰa̰ sʰɪ̀ɰ̃ dú]	ဒထွေး [da̰ dwé]	ဓအောက်ခြိုက် [da̰ ʔaʊʔ tɕʰaɪʔ]	နငယ် [na̰ ŋɛ̀]
Labials
(ဩဌဇ)
ပဝဂ်	ပ	p	/p/	ဖ	hp	/pʰ/	ဗ	b	/b/	ဘ	bh	/b/	မ	m	/m/
ပစောက် ([pa̰ zaʊʔ])	ဖဦးထုပ် ([pʰa̰ ʔóʊʔ tʰoʊʔ])	ဗထက်ခြိုက် ([ba̰ tɛʔ tɕʰaɪʔ])	ဘကုန်း ([ba̰ ɡóʊɰ̃])	မ [ma̰]
Miscellaneous consonants
Without group
(အဝဂ်)	ယ	y	/j/	ရ	r	/j/	လ	l	/l/	ဝ	w	/w/	သ	s	/θ/
ယပက်လက် [ja̰ pɛʔ lɛʔ]	ရကောက်‌ [ja̰ ɡaʊʔ]	လငယ် [la̰ ŋɛ̀]	ဝ‌ [wa̰]	သ‌ [θa̰]
ဟ	h	/h/	ဠ	l	/l/	အ	a	/ʔ/
ဟ‌ [ha̰]	ဠကြီး [la̰ dʑí]	အ [ʔa̰]
Independent vowels
ဣ	i.	/ʔḭ/	ဤ	i	/ʔì/	ဥ	u.	/ʔṵ/	ဦ	u	/ʔù/
ဧ	e	/ʔè/	ဩ	au:	/ʔɔ́/	ဪ
