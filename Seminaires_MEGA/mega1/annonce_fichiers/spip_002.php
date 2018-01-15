// Show or hide the element with the given 'id':
function refbase_toggleVisibilitySlide(id, imgid, txtid, txt, upd) {
	var e = document.getElementById(id);
	var i = document.getElementById(imgid);
	//var t = document.getElementById(txtid);
	//var upd;
	//if (upd === undefined) upd = true;
	if (e.style.display == 'block' || e.style.display == '') {
		//if (upd) 
		e.style.display = 'none';
		i.src = 'http://lmfa.ec-lyon.fr/plugins/refbase/refbase_closed.gif';
		//t.innerHTML = txt;
	}
	else {
		//if (upd)
		e.style.display = 'block';
		i.src = 'http://lmfa.ec-lyon.fr/plugins/refbase/refbase_open.gif';
		//t.innerHTML = '';
	}
}