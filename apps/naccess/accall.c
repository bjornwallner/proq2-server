/* accall.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__5000 = 5000;
static integer c__100 = 100;
static integer c__4 = 4;
static integer c__10 = 10;

/* accall_start */
/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static real hyrad = 1.f;
    static char rlab[3*3] = "RES" "HEM" "HOH";
    static integer nbackbone = 4;

    /* Format strings */
    static char fmt_102[] = "(\002 UNKNOWN residue type.............> \002,a"
	    "10)";
    static char fmt_104[] = "(\002 NON-STANDARD atom.|\002,a4,\002| in resid"
	    "ue> \002,a10)";
    static char fmt_106[] = "(\002 ASSUMED vdw of \002,a4,\002 in \002,a10"
	    ",\002 = \002,f5.2,\002 (same as \002,a3,\002)\002)";
    static char fmt_108[] = "(\002 GUESSED vdw of \002,a4,\002 in \002,a10"
	    ",\002 = \002,f5.2)";
    static char fmt_110[] = "(\002 CHAINS   \002,i5,/,\002 RESIDUES \002,i5,"
	    "/,\002 ATOMS    \002,i5)";
    static char fmt_120[] = "(\002REM  File of summed (Sum) and % (per."
	    ")\002,\002 accessibilities for \002,a)";
    static char fmt_125[] = "(\002REM RES _ NUM      All-atoms   Non-P-sid"
	    "e  \002,\002 Polar-Side   Total-Side   Main-Chain    \002,\002No"
	    "n-polar    All polar\002)";
    static char fmt_130[] = "(\002REM                ABS   REL    ABS   RE"
	    "L\002,\002    ABS   REL    ABS   REL    ABS   REL\002,\002    AB"
	    "S   REL    ABS   REL\002)";
    static char fmt_126[] = "(\002REM RES _ NUM      All-atoms   Total-Sid"
	    "e  \002,\002 Main-Chain    Non-polar    All polar\002)";
    static char fmt_131[] = "(\002REM                ABS   REL    ABS   RE"
	    "L\002,\002    ABS   REL    ABS   REL    ABS   REL\002)";
    static char fmt_150[] = "(a3,1x,a10,1x,7(f7.2,f6.1))";
    static char fmt_151[] = "(a3,1x,a10,1x,5(f7.2,f6.1))";
    static char fmt_154[] = "(\002END  Absolute sums over single chains surf"
	    "ace \002)";
    static char fmt_155[] = "(\002CHAIN \002,i2,1x,a1,3x,7(f8.1,5x))";
    static char fmt_156[] = "(\002CHAIN \002,i2,1x,a1,3x,5(f8.1,5x))";
    static char fmt_160[] = "(\002END  Absolute sums over all chains \002,/"
	    ",\002TOTAL\002,8x,7(f8.1,5x))";
    static char fmt_161[] = "(\002END  Absolute sums over all chains \002,/"
	    ",\002TOTAL\002,8x,5(f8.1,5x))";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
    real r__1, r__2;
    char ch__1[260];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_indx(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_stop(char *, ftnlen), s_cat(char *, char **, 
	    integer *, integer *, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), s_rsfi(icilist *), e_rsfi(void), f_clos(
	    cllist *);

    /* Local variables */
    static integer backbone[50000], resindex[50000];
    extern /* Subroutine */ int polguess_(char *, integer *, ftnlen);
    static integer polstats[50000];
    static char c__[256*256];
    static integer i__, j, k, l[256], n;
    extern doublereal readfloat_(char *, integer *, ftnlen);
    extern integer what_atom__(char *, integer *, integer *, ftnlen);
    static integer ip;
    static logical ok;
    static integer ir, num_chains__;
    extern integer readstring_(integer *, char *, integer *, ftnlen);
    static logical aok;
    static char alt[1], res[3];
    static real vdw;
    static integer rty[5000];
    static real xyz[150000]	/* was [50000][3] */, accs[50000];
    static char card[256];
    static integer flen, ilen;
    static real rads[50000];
    static logical falt;
    static char atom[4];
    static integer slen;
    static logical oldr;
    static char last[10];
    static integer vlen, nats;
    static real bfact[50000];
    static char label[30*50000];
    extern integer chain_(integer *, char *, char *, ftnlen, ftnlen);
    static char chnam[1*20], fname[256], sname[256];
    static logical hetas, conta;
    static char vname[256];
    extern integer fopen_(integer *, char *, integer *, char *, ftnlen, 
	    ftnlen);
    static real occup[50000];
    extern integer parse_(char *, integer *, char *, char *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static real probe;
    static logical dorsa;
    extern /* Subroutine */ int vanin_(char *, integer *, integer *, char *, 
	    char *, integer *, real *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), gatom_(char *, char *, integer *, integer *, logical *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical fullo;
    static integer atype;
    static logical resok;
    extern /* Subroutine */ int ratom_(char *, char *, integer *, integer *, 
	    logical *, integer *, ftnlen, ftnlen);
    static logical hydro;
    extern /* Subroutine */ int solva_(integer *, real *, real *, real *, 
	    real *, real *);
    static real csums[7];
    static logical start;
    static integer rtype[100];
    extern /* Subroutine */ int tolow_(char *, integer *, ftnlen);
    static real tsums[7];
    extern /* Subroutine */ int which3_(char *, char *, integer *, integer *, 
	    logical *, ftnlen, ftnlen);
    static integer achain[50000];
    static char aacids[3*100];
    static integer nchain, nacids, rchain[5000];
    static char anames[4*100*100];
    static real vradii[10000]	/* was [100][100] */;
    static char resnam[10*5000], firsta[1];
    static real zslice;
    static integer spolar[10000]	/* was [100][100] */, numats[100];
    extern /* Subroutine */ int summer_(char *, integer *, integer *, real *, 
	    integer *, integer *, integer *, integer *, char *, real *, real *
	    , ftnlen, ftnlen), vguess_(char *, real *, ftnlen);
    static integer num_res__;
    extern integer restype_(char *, ftnlen);
    static logical wwaters;
    static real ressums[70000]	/* was [5000][7][2] */;

    /* Fortran I/O blocks */
    static cilist io___35 = { 0, 4, 0, "(a)", 0 };
    static cilist io___36 = { 0, 4, 0, "(a,i7)", 0 };
    static cilist io___37 = { 0, 4, 0, "(a,i7)", 0 };
    static cilist io___38 = { 0, 4, 0, "(2a)", 0 };
    static cilist io___39 = { 0, 4, 0, "(a,f6.2)", 0 };
    static cilist io___40 = { 0, 4, 0, "(a,f6.3)", 0 };
    static cilist io___41 = { 0, 4, 0, "(2a)", 0 };
    static cilist io___42 = { 0, 4, 0, "(a)", 0 };
    static cilist io___43 = { 0, 4, 0, "(a)", 0 };
    static cilist io___44 = { 0, 4, 0, "(a)", 0 };
    static cilist io___45 = { 0, 4, 0, "(a)", 0 };
    static cilist io___46 = { 0, 4, 0, "(a)", 0 };
    static cilist io___47 = { 0, 4, 0, "(a)", 0 };
    static cilist io___48 = { 0, 4, 0, "(a,i3,a)", 0 };
    static cilist io___62 = { 0, 4, 0, "(a,i)", 0 };
    static cilist io___63 = { 0, 4, 0, "(a)", 0 };
    static cilist io___68 = { 0, 4, 0, fmt_102, 0 };
    static cilist io___75 = { 0, 4, 0, fmt_104, 0 };
    static cilist io___76 = { 0, 4, 0, fmt_106, 0 };
    static cilist io___77 = { 0, 4, 0, fmt_108, 0 };
    static icilist io___81 = { 0, card+54, 0, "(f)", 6, 1 };
    static icilist io___83 = { 0, card+60, 0, "(f)", 6, 1 };
    static icilist io___86 = { 0, card, 0, "(30x,3f8.3)", 256, 1 };
    static cilist io___88 = { 0, 4, 0, "(a)", 0 };
    static cilist io___89 = { 0, 4, 0, fmt_110, 0 };
    static cilist io___91 = { 0, 2, 0, "(a30,3f8.3,f6.2,f5.1,f7.3,1x,f5.2)", 
	    0 };
    static cilist io___92 = { 0, 2, 0, "(a30,3f8.3,f8.3,1x,f5.2)", 0 };
    static cilist io___93 = { 0, 4, 0, "(a)", 0 };
    static cilist io___96 = { 0, 4, 0, "(a)", 0 };
    static cilist io___97 = { 0, 3, 0, fmt_120, 0 };
    static cilist io___98 = { 0, 3, 0, fmt_125, 0 };
    static cilist io___99 = { 0, 3, 0, fmt_130, 0 };
    static cilist io___100 = { 0, 3, 0, fmt_126, 0 };
    static cilist io___101 = { 0, 3, 0, fmt_131, 0 };
    static cilist io___102 = { 0, 3, 0, fmt_150, 0 };
    static cilist io___103 = { 0, 3, 0, fmt_151, 0 };
    static cilist io___104 = { 0, 3, 0, fmt_154, 0 };
    static cilist io___106 = { 0, 3, 0, fmt_155, 0 };
    static cilist io___107 = { 0, 3, 0, fmt_156, 0 };
    static cilist io___108 = { 0, 3, 0, fmt_160, 0 };
    static cilist io___109 = { 0, 3, 0, fmt_161, 0 };



/*     --  VERSION: 2.1 */
/*     --  AIM:- */
/*     --  Input a Brookhaven entry file and output an */
/*     --  PDB format file after filtering/cleaning, including Van der */
/*     --  Waal radii, contained in an external file "vdw.radii". */

/*     --  INPUT:- */
/*     --  PDB format file, van der Waal radii file */

/*     --  OPTIONS:- */
/*     --  Inclusion of non-standard amino acids, het-atoms, waters */
/*     --  etc. Flagging of missing residues, chain-breaks, */
/*     --  non-standard atoms names, missing atoms. */
/*     --  nucleic acids, separable chains, polar/non-polar summing */

/*     --  AUTHOR: S. Hubbard 3/92. EMBL. */


/*     -- functions */


/* Parameters uses in ACCALL.F */
/*     --  PARAMETERS: maxa = max atoms per residue */
/*     --	      maxr = max number of "standard" residue types */
/*     --              maxs = max number of atoms in PDB file */
/*     --              maxx = max number of residues in PDB file */
/*     --              maxc = max number of chains */

/* ============================================================= */

/* SOLVA LOCAL PARAMETERS: */
/* ncube = maximum number of cubes allowed for placing of atoms */
/* nac   = maximum number of atoms per cube */
/* nint  = maximum number of sphere intersections */


/*     -- variables */


/*     -- defaults */

    hetas = FALSE_;
    hydro = FALSE_;
    wwaters = FALSE_;
    fullo = FALSE_;
    dorsa = TRUE_;
    conta = FALSE_;
    oldr = FALSE_;

/*     -- Get USER directives */

    while(readstring_(&c__5, card, &ilen, (ftnlen)256) >= 0) {
	n = parse_(card, &ilen, " ", c__, l, (ftnlen)256, (ftnlen)1, (ftnlen)
		256);
	tolow_(c__, l, (ftnlen)256);
	if (s_cmp(c__, "pdbf", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(fname, c__ + 256, (ftnlen)256, (ftnlen)256);
	    flen = l[1];
	} else if (s_cmp(c__, "vdwf", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(vname, c__ + 256, (ftnlen)256, (ftnlen)256);
	    vlen = l[1];
	} else if (s_cmp(c__, "stdf", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(sname, c__ + 256, (ftnlen)256, (ftnlen)256);
	    slen = l[1];
	} else if (s_cmp(c__, "prob", (ftnlen)4, (ftnlen)4) == 0) {
	    probe = readfloat_(c__ + 256, &l[1], (ftnlen)256);
	} else if (s_cmp(c__, "zsli", (ftnlen)4, (ftnlen)4) == 0) {
	    zslice = readfloat_(c__ + 256, &l[1], (ftnlen)256);
	} else if (s_cmp(c__, "heta", (ftnlen)4, (ftnlen)4) == 0) {
	    hetas = TRUE_;
	} else if (s_cmp(c__, "hydr", (ftnlen)4, (ftnlen)4) == 0) {
	    hydro = TRUE_;
	} else if (s_cmp(c__, "wate", (ftnlen)4, (ftnlen)4) == 0) {
	    wwaters = TRUE_;
	} else if (s_cmp(c__, "full", (ftnlen)4, (ftnlen)4) == 0) {
	    fullo = TRUE_;
	} else if (s_cmp(c__, "oldr", (ftnlen)4, (ftnlen)4) == 0) {
	    oldr = TRUE_;
	} else if (s_cmp(c__, "asao", (ftnlen)4, (ftnlen)4) == 0) {
	    dorsa = FALSE_;
	} else if (s_cmp(c__, "cont", (ftnlen)4, (ftnlen)4) == 0) {
	    conta = TRUE_;
	} else if (s_cmp(c__, "csid", (ftnlen)4, (ftnlen)4) == 0) {
	    nbackbone = 5;
	}
    }

/*     --   open files */

    i__ = i_indx(fname, ".", (ftnlen)256, (ftnlen)1) - 1;
    k = 1;
    ok = FALSE_;
    for (j = i__; j >= 1; --j) {
	if (*(unsigned char *)&fname[j - 1] == '/' && ! ok) {
	    k = j + 1;
	    ok = TRUE_;
	}
    }
    if (fopen_(&c__1, fname, &flen, "old", (ftnlen)256, (ftnlen)3) == 0) {
	s_stop("ERROR: unable to open PDB file", (ftnlen)30);
    }
    o__1.oerr = 0;
    o__1.ounit = 2;
    o__1.ofnmlen = i__ - (k - 1) + 4;
/* Writing concatenation */
    i__1[0] = i__ - (k - 1), a__1[0] = fname + (k - 1);
    i__1[1] = 4, a__1[1] = ".asa";
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)260);
    o__1.ofnm = ch__1;
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    if (dorsa) {
	o__1.oerr = 0;
	o__1.ounit = 3;
	o__1.ofnmlen = i__ - (k - 1) + 4;
/* Writing concatenation */
	i__1[0] = i__ - (k - 1), a__1[0] = fname + (k - 1);
	i__1[1] = 4, a__1[1] = ".rsa";
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)260);
	o__1.ofnm = ch__1;
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = i__ - (k - 1) + 4;
/* Writing concatenation */
    i__1[0] = i__ - (k - 1), a__1[0] = fname + (k - 1);
    i__1[1] = 4, a__1[1] = ".log";
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)260);
    o__1.ofnm = ch__1;
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/*     -- Read in VDW radii for all residues/atoms */

    vanin_(vname, &vlen, &nacids, aacids, anames, numats, vradii, spolar, 
	    rtype, (ftnlen)256, (ftnlen)3, (ftnlen)4);
    s_wsfe(&io___35);
    do_fio(&c__1, " ACCALL - Accessibility calculations", (ftnlen)36);
    e_wsfe();
    s_wsfe(&io___36);
    do_fio(&c__1, " MAX RESIDUES  ", (ftnlen)15);
    do_fio(&c__1, (char *)&c__5000, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___37);
    do_fio(&c__1, " MAX ATOMS/RES ", (ftnlen)15);
    do_fio(&c__1, (char *)&c__100, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___38);
    do_fio(&c__1, " PDB FILE INPUT ", (ftnlen)16);
    do_fio(&c__1, fname, flen);
    e_wsfe();
    s_wsfe(&io___39);
    do_fio(&c__1, " PROBE SIZE     ", (ftnlen)16);
    do_fio(&c__1, (char *)&probe, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___40);
    do_fio(&c__1, " Z-SLICE WIDTH  ", (ftnlen)16);
    do_fio(&c__1, (char *)&zslice, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___41);
    do_fio(&c__1, " VDW RADII FILE ", (ftnlen)16);
    do_fio(&c__1, vname, vlen);
    e_wsfe();
    if (hetas) {
	s_wsfe(&io___42);
	do_fio(&c__1, " INCL HETATOMS", (ftnlen)14);
	e_wsfe();
    } else {
	s_wsfe(&io___43);
	do_fio(&c__1, " EXCL HETATOMS", (ftnlen)14);
	e_wsfe();
    }
    if (hydro) {
	s_wsfe(&io___44);
	do_fio(&c__1, " INCL HYDROGENS", (ftnlen)15);
	e_wsfe();
    } else {
	s_wsfe(&io___45);
	do_fio(&c__1, " EXCL HYDROGENS", (ftnlen)15);
	e_wsfe();
    }
    if (wwaters) {
	s_wsfe(&io___46);
	do_fio(&c__1, " INCL WATERS", (ftnlen)12);
	e_wsfe();
    } else {
	s_wsfe(&io___47);
	do_fio(&c__1, " EXCL WATERS", (ftnlen)12);
	e_wsfe();
    }
    s_wsfe(&io___48);
    do_fio(&c__1, " READVDW ", (ftnlen)9);
    do_fio(&c__1, (char *)&nacids, (ftnlen)sizeof(integer));
    do_fio(&c__1, " residues input", (ftnlen)15);
    e_wsfe();

/* --  Initialise variables/logicals */

    falt = FALSE_;
    start = TRUE_;
    s_copy(last, "          ", (ftnlen)10, (ftnlen)10);

/* -- Read data & decode */

    nats = 0;
    num_res__ = 0;
    nchain = 1;
    while(readstring_(&c__1, card, &ilen, (ftnlen)256) >= 0) {
	atype = restype_(card, (ftnlen)20);
	ip = -1;
	if (atype == 1 || atype == 2 && hetas || atype == 3 && wwaters) {

/*     -- Ignore Alternate positions, other than blanks or 1st */
/*     -- encountered */

	    *(unsigned char *)alt = *(unsigned char *)&card[16];
	    if (*(unsigned char *)alt != ' ') {
		if (! falt) {
		    *(unsigned char *)firsta = *(unsigned char *)alt;
		    falt = TRUE_;
		}
		if (*(unsigned char *)alt != *(unsigned char *)firsta) {
		    goto L5;
		}
	    }

/*     -- Ignore hydrogens & deuteriums (unless flagged) */

	    if (*(unsigned char *)&card[13] == 'H' || *(unsigned char *)&card[
		    13] == 'D') {
		if (! hydro) {
		    goto L5;
		}
		vdw = hyrad;
		++nats;
		goto L6;
	    }

/*     -- Next atom */

	    ++nats;

/*     -- First residue ? */

	    if (start) {
		start = FALSE_;
		*(unsigned char *)&chnam[0] = *(unsigned char *)&card[21];
		num_chains__ = 1;
	    }

/*     -- New residue ? */

	    if (s_cmp(last, card + 17, (ftnlen)10, (ftnlen)10) != 0) {
		s_copy(last, card + 17, (ftnlen)10, (ftnlen)10);
		++num_res__;
		if (num_res__ > 5000) {
		    s_wsfe(&io___62);
		    do_fio(&c__1, " ERROR - Maximum number of residues excee"
			    "ded ", (ftnlen)45);
		    do_fio(&c__1, (char *)&c__5000, (ftnlen)sizeof(integer));
		    e_wsfe();
		    s_wsfe(&io___63);
		    do_fio(&c__1, " Increase maxx in accall.pars and recompi"
			    "le", (ftnlen)43);
		    e_wsfe();
		    s_stop("SOLVA_ERROR: maxx exceeded", (ftnlen)26);
		}
		s_copy(res, card + 17, (ftnlen)3, (ftnlen)3);
		s_copy(resnam + (num_res__ - 1) * 10, last, (ftnlen)10, (
			ftnlen)10);
		which3_(res, aacids, &i__, &nacids, &resok, (ftnlen)3, (
			ftnlen)3);
		if (i__ > 0) {
		    ir = rtype[i__ - 1];
		} else {
		    ir = 0;
		}
		if (! resok) {
		    s_wsfe(&io___68);
		    do_fio(&c__1, card + 17, (ftnlen)10);
		    e_wsfe();
		}
		rty[num_res__ - 1] = atype;
		nchain = chain_(&num_chains__, chnam, card + 21, (ftnlen)1, (
			ftnlen)1);
		rchain[num_res__ - 1] = nchain;
	    }

/*     -- Get atom type */

	    s_copy(atom, card + 12, (ftnlen)4, (ftnlen)4);
	    backbone[nats - 1] = what_atom__(atom, &ir, &nbackbone, (ftnlen)4)
		    ;
	    achain[nats - 1] = nchain;

/*     -- Assign radius to atom */
/*     -- Special case(s) */

	    if (s_cmp(atom, " OXT", (ftnlen)4, (ftnlen)4) == 0) {
		vdw = 1.4f;
		ip = 1;
		goto L6;
	    }

/*     -- known residue type */

	    if (resok) {
		vdw = 0.f;
		ratom_(atom, anames, &i__, numats, &aok, &j, (ftnlen)4, (
			ftnlen)4);
	    }

/*     -- Not OK, then try atoms in all residue types */

	    if (! resok || ! aok) {
		gatom_(atom, anames, &nacids, numats, &aok, &i__, &j, (ftnlen)
			4, (ftnlen)4);
		s_wsfe(&io___75);
		do_fio(&c__1, atom, (ftnlen)4);
		do_fio(&c__1, card + 17, (ftnlen)10);
		e_wsfe();
		if (aok) {
		    s_wsfe(&io___76);
		    do_fio(&c__1, atom, (ftnlen)4);
		    do_fio(&c__1, card + 17, (ftnlen)10);
		    do_fio(&c__1, (char *)&vradii[i__ + j * 100 - 101], (
			    ftnlen)sizeof(real));
		    do_fio(&c__1, aacids + (i__ - 1) * 3, (ftnlen)3);
		    e_wsfe();
		}
	    }

/*     -- Still not OK, make a guess */

	    if (! aok) {
		vguess_(atom, &vdw, (ftnlen)4);
		s_wsfe(&io___77);
		do_fio(&c__1, atom, (ftnlen)4);
		do_fio(&c__1, card + 17, (ftnlen)10);
		do_fio(&c__1, (char *)&vdw, (ftnlen)sizeof(real));
		e_wsfe();
	    } else {
		vdw = vradii[i__ + j * 100 - 101];
		ip = spolar[i__ + j * 100 - 101];
	    }
	    if (ip < 0) {
		polguess_(atom, &ip, (ftnlen)4);
	    }

/*     -- Store data */

L6:
	    rads[nats - 1] = vdw;
	    polstats[nats - 1] = ip;
	    s_copy(label + (nats - 1) * 30, card, (ftnlen)30, (ftnlen)30);
	    if (fullo) {
		s_rsfi(&io___81);
		do_fio(&c__1, (char *)&occup[nats - 1], (ftnlen)sizeof(real));
		e_rsfi();
		s_rsfi(&io___83);
		do_fio(&c__1, (char *)&bfact[nats - 1], (ftnlen)sizeof(real));
		e_rsfi();
	    }
	    resindex[nats - 1] = num_res__;
	    s_rsfi(&io___86);
	    for (k = 1; k <= 3; ++k) {
		do_fio(&c__1, (char *)&xyz[nats + k * 50000 - 50001], (ftnlen)
			sizeof(real));
	    }
	    e_rsfi();
	}
L5:
	;
    }

/*     -- output */

    cl__1.cerr = 0;
    cl__1.cunit = 1;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_wsfe(&io___88);
    do_fio(&c__1, " ADDED VDW RADII", (ftnlen)16);
    e_wsfe();
    s_wsfe(&io___89);
    do_fio(&c__1, (char *)&num_chains__, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&num_res__, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nats, (ftnlen)sizeof(integer));
    e_wsfe();

/*     -- calculate atomic accessibilities */

    solva_(&nats, xyz, rads, accs, &probe, &zslice);
    if (conta) {
	i__2 = nats;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    r__1 = rads[i__ - 1];
/* Computing 2nd power */
	    r__2 = rads[i__ - 1] + probe;
	    accs[i__ - 1] = accs[i__ - 1] * (r__1 * r__1) / (r__2 * r__2);
	}
    }
    i__2 = nats;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (fullo) {
	    s_wsfe(&io___91);
	    do_fio(&c__1, label + (i__ - 1) * 30, (ftnlen)30);
	    for (j = 1; j <= 3; ++j) {
		do_fio(&c__1, (char *)&xyz[i__ + j * 50000 - 50001], (ftnlen)
			sizeof(real));
	    }
	    do_fio(&c__1, (char *)&occup[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&bfact[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&accs[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&rads[i__ - 1], (ftnlen)sizeof(real));
	    e_wsfe();
	} else {
	    s_wsfe(&io___92);
	    do_fio(&c__1, label + (i__ - 1) * 30, (ftnlen)30);
	    for (j = 1; j <= 3; ++j) {
		do_fio(&c__1, (char *)&xyz[i__ + j * 50000 - 50001], (ftnlen)
			sizeof(real));
	    }
	    do_fio(&c__1, (char *)&accs[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&rads[i__ - 1], (ftnlen)sizeof(real));
	    e_wsfe();
	}
    }
    s_wsfe(&io___93);
    do_fio(&c__1, " CALCULATED ATOMIC ACCESSIBILITES", (ftnlen)33);
    e_wsfe();
    if (dorsa) {
	summer_(sname, &slen, &nats, accs, backbone, polstats, rtype, 
		resindex, resnam, ressums, tsums, (ftnlen)256, (ftnlen)10);
	s_wsfe(&io___96);
	do_fio(&c__1, " SUMMED ACCESSIBILITIES OVER RESIDUES", (ftnlen)37);
	e_wsfe();
	s_wsfe(&io___97);
	e_wsfe();
	if (oldr) {
	    s_wsfe(&io___98);
	    e_wsfe();
	    s_wsfe(&io___99);
	    e_wsfe();
	} else {
	    s_wsfe(&io___100);
	    e_wsfe();
	    s_wsfe(&io___101);
	    e_wsfe();
	}
	i__2 = resindex[nats - 1];
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (oldr) {
		s_wsfe(&io___102);
		do_fio(&c__1, rlab + (rty[i__ - 1] - 1) * 3, (ftnlen)3);
		do_fio(&c__1, resnam + (i__ - 1) * 10, (ftnlen)10);
		for (j = 1; j <= 7; ++j) {
		    do_fio(&c__1, (char *)&ressums[i__ + (j + 7) * 5000 - 
			    40001], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&ressums[i__ + (j + 14) * 5000 - 
			    40001], (ftnlen)sizeof(real));
		}
		e_wsfe();
	    } else {
		s_wsfe(&io___103);
		do_fio(&c__1, rlab + (rty[i__ - 1] - 1) * 3, (ftnlen)3);
		do_fio(&c__1, resnam + (i__ - 1) * 10, (ftnlen)10);
		do_fio(&c__1, (char *)&ressums[i__ - 1], (ftnlen)sizeof(real))
			;
		do_fio(&c__1, (char *)&ressums[i__ + 34999], (ftnlen)sizeof(
			real));
		for (j = 4; j <= 7; ++j) {
		    do_fio(&c__1, (char *)&ressums[i__ + (j + 7) * 5000 - 
			    40001], (ftnlen)sizeof(real));
		    do_fio(&c__1, (char *)&ressums[i__ + (j + 14) * 5000 - 
			    40001], (ftnlen)sizeof(real));
		}
		e_wsfe();
	    }
	}
	s_wsfe(&io___104);
	e_wsfe();
	i__2 = num_chains__;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (*(unsigned char *)&chnam[i__ - 1] == ' ') {
		*(unsigned char *)&chnam[i__ - 1] = '_';
	    }
	    for (j = 1; j <= 7; ++j) {
		csums[j - 1] = 0.f;
	    }
	    i__3 = num_res__;
	    for (j = 1; j <= i__3; ++j) {
		if (rchain[j - 1] == i__) {
		    for (k = 1; k <= 7; ++k) {
			csums[k - 1] += ressums[j + (k + 7) * 5000 - 40001];
		    }
		}
	    }
	    if (oldr) {
		s_wsfe(&io___106);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, chnam + (i__ - 1), (ftnlen)1);
		for (j = 1; j <= 7; ++j) {
		    do_fio(&c__1, (char *)&csums[j - 1], (ftnlen)sizeof(real))
			    ;
		}
		e_wsfe();
	    } else {
		s_wsfe(&io___107);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_fio(&c__1, chnam + (i__ - 1), (ftnlen)1);
		do_fio(&c__1, (char *)&csums[0], (ftnlen)sizeof(real));
		for (j = 4; j <= 7; ++j) {
		    do_fio(&c__1, (char *)&csums[j - 1], (ftnlen)sizeof(real))
			    ;
		}
		e_wsfe();
	    }
	}
	if (oldr) {
	    s_wsfe(&io___108);
	    for (i__ = 1; i__ <= 7; ++i__) {
		do_fio(&c__1, (char *)&tsums[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	} else {
	    s_wsfe(&io___109);
	    do_fio(&c__1, (char *)&tsums[0], (ftnlen)sizeof(real));
	    for (i__ = 4; i__ <= 7; ++i__) {
		do_fio(&c__1, (char *)&tsums[i__ - 1], (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	}
    }
/* ---------------------------------------------------------- */
/* ---------------------------------------------------------- */
    return 0;
} /* MAIN__ */

/* Subroutine */ int vguess_(char *atom, real *vdw, ftnlen atom_len)
{
    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    *vdw = 1.8f;

/* -- Make a guess then ! */

    if (*(unsigned char *)&atom[1] == 'C') {
	*vdw = 1.8f;
    }
    if (*(unsigned char *)&atom[1] == 'N') {
	*vdw = 1.6f;
    }
    if (*(unsigned char *)&atom[1] == 'S') {
	*vdw = 1.85f;
    }
    if (*(unsigned char *)&atom[1] == 'O') {
	*vdw = 1.4f;
    }
    if (*(unsigned char *)&atom[1] == 'P') {
	*vdw = 1.9f;
    }
    if (s_cmp(atom, "CA", (ftnlen)2, (ftnlen)2) == 0) {
	*vdw = 2.07f;
    }
    if (s_cmp(atom, "FE", (ftnlen)2, (ftnlen)2) == 0) {
	*vdw = 1.47f;
    }
    if (s_cmp(atom, "CU", (ftnlen)2, (ftnlen)2) == 0) {
	*vdw = 1.78f;
    }
    if (s_cmp(atom, "ZN", (ftnlen)2, (ftnlen)2) == 0) {
	*vdw = 1.39f;
    }
    if (s_cmp(atom, "MG", (ftnlen)2, (ftnlen)2) == 0) {
	*vdw = 1.73f;
    }
    return 0;
} /* vguess_ */


/*     Return integer number of chain, assigned from */
/*     the single letter. This copes with different */
/*     parts of the chain in different parts of the PDB file. */

integer chain_(integer *nc, char *names, char *c__, ftnlen names_len, ftnlen 
	c_len)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__;


/* Parameters uses in ACCALL.F */
/*     --  PARAMETERS: maxa = max atoms per residue */
/*     --	      maxr = max number of "standard" residue types */
/*     --              maxs = max number of atoms in PDB file */
/*     --              maxx = max number of residues in PDB file */
/*     --              maxc = max number of chains */

/* ============================================================= */

/* SOLVA LOCAL PARAMETERS: */
/* ncube = maximum number of cubes allowed for placing of atoms */
/* nac   = maximum number of atoms per cube */
/* nint  = maximum number of sphere intersections */

    /* Parameter adjustments */
    --names;

    /* Function Body */
    i__1 = *nc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)c__ == *(unsigned char *)&names[i__]) {
	    ret_val = i__;
	    return ret_val;
	}
    }
    ++(*nc);
    *(unsigned char *)&names[*nc] = *(unsigned char *)c__;
    ret_val = *nc;
    return ret_val;
} /* chain_ */


/*     Subroutine Gatom finds the residue and atom index */
/*     (the latter from ratom, see below) */


/* Subroutine */ int gatom_(char *atom, char *anames, integer *nres, integer *
	nats, logical *ok, integer *ir, integer *ia, ftnlen atom_len, ftnlen 
	anames_len)
{
    extern /* Subroutine */ int ratom_(char *, char *, integer *, integer *, 
	    logical *, integer *, ftnlen, ftnlen);


/* Parameters uses in ACCALL.F */
/*     --  PARAMETERS: maxa = max atoms per residue */
/*     --	      maxr = max number of "standard" residue types */
/*     --              maxs = max number of atoms in PDB file */
/*     --              maxx = max number of residues in PDB file */
/*     --              maxc = max number of chains */

/* ============================================================= */

/* SOLVA LOCAL PARAMETERS: */
/* ncube = maximum number of cubes allowed for placing of atoms */
/* nac   = maximum number of atoms per cube */
/* nint  = maximum number of sphere intersections */

    /* Parameter adjustments */
    --nats;
    anames -= 404;

    /* Function Body */
    *ok = FALSE_;
    *ir = 0;
    while(*ir < *nres && ! (*ok)) {
	++(*ir);
	ratom_(atom, anames + 404, ir, &nats[1], ok, ia, (ftnlen)4, (ftnlen)4)
		;
    }
    return 0;
} /* gatom_ */

/* Subroutine */ int ratom_(char *atom, char *anames, integer *ires, integer *
	nats, logical *ok, integer *find, ftnlen atom_len, ftnlen anames_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;


/* --  Checks to see if a "standard" atom name has been read in. */
/* --  Standard atom names are read from the file "vdw.radii", for */
/* --  the defined residue types therein. OK=.true. if found, and */
/* --  residue type = ires, and find = integer identifier of atom. */


/* Parameters uses in ACCALL.F */
/*     --  PARAMETERS: maxa = max atoms per residue */
/*     --	      maxr = max number of "standard" residue types */
/*     --              maxs = max number of atoms in PDB file */
/*     --              maxx = max number of residues in PDB file */
/*     --              maxc = max number of chains */

/* ============================================================= */

/* SOLVA LOCAL PARAMETERS: */
/* ncube = maximum number of cubes allowed for placing of atoms */
/* nac   = maximum number of atoms per cube */
/* nint  = maximum number of sphere intersections */

    /* Parameter adjustments */
    --nats;
    anames -= 404;

    /* Function Body */
    *ok = FALSE_;
    *find = 0;
    if (*ires == 0 || *ires > 100) {
	return 0;
    }
    i__1 = nats[*ires];
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s_cmp(atom, anames + (*ires + i__ * 100 << 2), (ftnlen)4, (ftnlen)
		4) == 0) {
	    *find = i__;
	    *ok = TRUE_;
	    return 0;
	}
    }
    return 0;
} /* ratom_ */

integer restype_(char *card, ftnlen card_len)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    ret_val = 0;
    if (s_cmp(card, "ATOM", (ftnlen)4, (ftnlen)4) == 0) {
	ret_val = 1;
    }
    if (s_cmp(card, "HETATM", (ftnlen)6, (ftnlen)6) == 0) {
	ret_val = 2;
    }
    if (s_cmp(card + 17, "HOH", (ftnlen)3, (ftnlen)3) == 0) {
	ret_val = 3;
    }
    return ret_val;
} /* restype_ */

/* Subroutine */ int which3_(char *res, char *acids, integer *ires, integer *
	nacids, logical *ok, ftnlen res_len, ftnlen acids_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;


/* -- Search array "acids" for existence of residue "res". */
/* -- OK = .true. if found, and index of res returned in "ires". */


/* Parameters uses in ACCALL.F */
/*     --  PARAMETERS: maxa = max atoms per residue */
/*     --	      maxr = max number of "standard" residue types */
/*     --              maxs = max number of atoms in PDB file */
/*     --              maxx = max number of residues in PDB file */
/*     --              maxc = max number of chains */

/* ============================================================= */

/* SOLVA LOCAL PARAMETERS: */
/* ncube = maximum number of cubes allowed for placing of atoms */
/* nac   = maximum number of atoms per cube */
/* nint  = maximum number of sphere intersections */

    /* Parameter adjustments */
    acids -= 3;

    /* Function Body */
    *ires = *nacids + 1;
    i__1 = *nacids;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s_cmp(res, acids + i__ * 3, (ftnlen)3, (ftnlen)3) == 0) {
	    *ires = i__;
	    *ok = TRUE_;
	    return 0;
	}
    }
    *ok = FALSE_;
    return 0;
} /* which3_ */


/*     Open a file to i/o channel "iochan" */

integer fopen_(integer *iochan, char *filename, integer *flen, char *fstat, 
	ftnlen filename_len, ftnlen fstat_len)
{
    /* System generated locals */
    integer ret_val, i__1;
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *);

    o__1.oerr = 1;
    o__1.ounit = *iochan;
    o__1.ofnmlen = *flen;
    o__1.ofnm = filename;
    o__1.orl = 0;
    o__1.osta = fstat;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L100;
    }
    ret_val = 1;
    return ret_val;
L100:
    ret_val = 0;
    return ret_val;
} /* fopen_ */

integer readstring_(integer *file, char *card, integer *flen, ftnlen card_len)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void);

    /* Fortran I/O blocks */
    static cilist io___113 = { 1, 0, 1, "(a)", 0 };


    if (*file > 200) {
	s_stop("ERROR: file number too large", (ftnlen)28);
    }
    io___113.ciunit = *file;
    i__1 = s_rsfe(&io___113);
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = do_fio(&c__1, card, (ftnlen)256);
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L100;
    }
    *flen = 256;
    while(*(unsigned char *)&card[*flen - 1] == ' ') {
	--(*flen);
    }
    ret_val = *flen;
    return ret_val;
L100:
    ret_val = -1;
    return ret_val;
} /* readstring_ */

integer parse_(char *card, integer *length, char *separator, char *chars, 
	integer *clen, ftnlen card_len, ftnlen separator_len, ftnlen 
	chars_len)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static logical search;
    static integer charpos[512]	/* was [256][2] */;

    /* Parameter adjustments */
    --clen;
    chars -= chars_len;

    /* Function Body */
    ret_val = 0;
    i__ = 0;
    search = FALSE_;
    while(i__ < *length) {
	++i__;
	if (! search) {
	    if (*(unsigned char *)&card[i__ - 1] != *(unsigned char *)
		    separator) {
		++ret_val;
		charpos[ret_val - 1] = i__;
		search = TRUE_;
	    }
	} else {
	    if (*(unsigned char *)&card[i__ - 1] == *(unsigned char *)
		    separator) {
		charpos[ret_val + 255] = i__ - 1;
		search = FALSE_;
	    }
	}
    }
    if (search) {
	charpos[ret_val + 255] = *length;
    }
    i__1 = ret_val;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = charpos[i__ - 1] - 1;
	s_copy(chars + i__ * chars_len, card + i__2, chars_len, charpos[i__ + 
		255] - i__2);
	clen[i__] = charpos[i__ + 255] - charpos[i__ - 1] + 1;
    }
    return ret_val;
} /* parse_ */

integer readint_(char *card, integer *clen, ftnlen card_len)
{
    /* System generated locals */
    integer ret_val, i__1;
    icilist ici__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void);

    /* Local variables */
    static real value;

    ici__1.icierr = 1;
    ici__1.iciend = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = *clen;
    ici__1.iciunit = card;
    ici__1.icifmt = 0;
    i__1 = s_rsli(&ici__1);
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&value, (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = e_rsli();
    if (i__1 != 0) {
	goto L100;
    }
    ret_val = (integer) value;
    return ret_val;
L100:
    ret_val = -999;
    return ret_val;
} /* readint_ */

doublereal readfloat_(char *card, integer *clen, ftnlen card_len)
{
    /* System generated locals */
    integer i__1;
    real ret_val;
    icilist ici__1;

    /* Builtin functions */
    integer s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void);

    /* Local variables */
    static real value;

    ici__1.icierr = 1;
    ici__1.iciend = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = *clen;
    ici__1.iciunit = card;
    ici__1.icifmt = 0;
    i__1 = s_rsli(&ici__1);
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&value, (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = e_rsli();
    if (i__1 != 0) {
	goto L100;
    }
    ret_val = value;
    return ret_val;
L100:
    ret_val = -999.9f;
    return ret_val;
} /* readfloat_ */


/*     Convert uppercase to lowercase characters */

/* Subroutine */ int tolow_(char *text, integer *clen, ftnlen text_len)
{
    /* Initialized data */

    static char alph[1*26*2] = "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" 
	    "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y" "Z" "a" 
	    "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "m" "n" "o" "p" "q" 
	    "r" "s" "t" "u" "v" "w" "x" "y" "z";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static logical ok;

    i__1 = *clen;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = 0;
	ok = TRUE_;
	while(j < 26 && ok) {
	    ++j;
	    if (*(unsigned char *)&alph[j - 1] == *(unsigned char *)&text[i__ 
		    - 1]) {
		ok = FALSE_;
		*(unsigned char *)&text[i__ - 1] = *(unsigned char *)&alph[j 
			+ 25];
	    }
	}
    }
    return 0;
} /* tolow_ */

/* Subroutine */ int solva_(integer *nats, real *xyz, real *rads, real *accs, 
	real *probe, real *zslice)
{
    /* Initialized data */

    static real xmin = 9999.f;
    static real ymin = 9999.f;
    static real zmin = 9999.f;
    static real xmax = -9999.f;
    static real ymax = -9999.f;
    static real zmax = -9999.f;

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2;

    /* Builtin functions */
    double acos(doublereal);
    /* Subroutine */ int s_stop(char *, ftnlen);
    double sqrt(doublereal), atan2(doublereal, doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static real b, d__[2000];
    static integer i__, j, k, l, m, n;
    static real t, trig_test__;
    static integer in, io;
    static real pi, tf;
    static integer ir;
    static real dx[2000], dy[2000];
    static integer nm;
    static real ti, rr, tt, xr, yr, zr, rad[50000];
    static integer tag[2000], kji, ict;
    static real dsq[2000];
    static integer nzp;
    static real pix2, rrx2, area, arcf[2000], beta;
    static integer cube[50000];
    static real arci[2000];
    static integer itab[9000000], karc, idim, mkji, natm[1350000000]	/* 
	    was [150][9000000] */;
    static real rmax;
    static integer inov[2000];
    static real zres, rrsq, alpha, parea;
    static integer jidim;
    static real radsq[50000], rsecn, rsecr, zgrid, rsec2n, rsec2r;
    static integer kjidim;
    static real arcsum;
    extern /* Subroutine */ int sortag_(real *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___186 = { 0, 4, 0, "(a)", 0 };


/* ** */
/* ************************************************************* */
/* **  SOLVA - LEE & RICHARDS TYPE ACCESSIBLITY CALCULATIONS ** */
/* ************************************************************* */
/* ** */
/* ** Calculate accessible surface area for a group of atoms. */
/* ** The accessible area for a given atom is calculated by the */
/* ** formula: */
/* **     (arcsum) x (atom radius+probe radius) x (deltaz) */
/* ** Numerical integration is carried out over z. in each z- */
/* ** section, the arcsum for a given atom is the arclength of */
/* ** the circle (intersection of the atom sphere with the z- */
/* ** section) that is not interior to any other atom circles */
/* ** in the same z-section. */
/* ** */
/* ************************************************************* */
/* ** */
/* **  error parameter  - this gives accuracy of calculation */
/* **                   - suitable values are 0.01 (high */
/* **                     accuracy) to 0.1 (low accuracy) */
/* **                   - in detail the z sections are spaced */
/* **                     at about error*diameter of atom */
/* ** */
/* **  probe size       - radius of probe in angstroms */
/* **                   - suitable value for water = 1.4 */
/* ** */
/* ************************************************************* */

/*     the following are dimensioned to the max no of atoms (maxs) */


/* Parameters uses in ACCALL.F */
/*     --  PARAMETERS: maxa = max atoms per residue */
/*     --	      maxr = max number of "standard" residue types */
/*     --              maxs = max number of atoms in PDB file */
/*     --              maxx = max number of residues in PDB file */
/*     --              maxc = max number of chains */

/* ============================================================= */

/* SOLVA LOCAL PARAMETERS: */
/* ncube = maximum number of cubes allowed for placing of atoms */
/* nac   = maximum number of atoms per cube */
/* nint  = maximum number of sphere intersections */


/*     the following are dimensioned to the max no of intersections */
/*     of neighbouring spheres (nint) */

    /* Parameter adjustments */
    --accs;
    --rads;
    xyz -= 50001;

    /* Function Body */

/*     initialise variables, constants */

    ict = 2000;
    pi = acos(-1.f);
    pix2 = pi * 2.f;

/*     -- Radius of an atom sphere = atom radius + probe radius */
/*     -- Find maxima and minima */

    rmax = 0.f;
    karc = ict;
    i__1 = *nats;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rad[i__ - 1] = rads[i__] + *probe;
/* Computing 2nd power */
	r__1 = rad[i__ - 1];
	radsq[i__ - 1] = r__1 * r__1;
	if (rad[i__ - 1] > rmax) {
	    rmax = rad[i__ - 1];
	}
	if (xmin > xyz[i__ + 50000]) {
	    xmin = xyz[i__ + 50000];
	}
	if (ymin > xyz[i__ + 100000]) {
	    ymin = xyz[i__ + 100000];
	}
	if (zmin > xyz[i__ + 150000]) {
	    zmin = xyz[i__ + 150000];
	}
	if (xmax < xyz[i__ + 50000]) {
	    xmax = xyz[i__ + 50000];
	}
	if (ymax < xyz[i__ + 100000]) {
	    ymax = xyz[i__ + 100000];
	}
	if (zmax < xyz[i__ + 150000]) {
	    zmax = xyz[i__ + 150000];
	}
    }

/*     rmax = max diameter */

    rmax *= 2.f;

/*     -- Cubicals containing the atoms are setup. */
/*     -- The dimension of an edge equals the largest atom sphere radius */
/*     -- The cubes have a single index */
/*     -- Minimum of 3 by 3 cubic grid */
/*     -- EXIT if max cubes exceeded */

    idim = (xmax - xmin) / rmax + 1.f;
    if (idim < 3) {
	idim = 3;
    }
    jidim = (ymax - ymin) / rmax + 1.f;
    if (jidim < 3) {
	jidim = 3;
    }
    jidim = idim * jidim;
    kjidim = (zmax - zmin) / rmax + 1.f;
    if (kjidim < 3) {
	kjidim = 3;
    }
    kjidim = jidim * kjidim;
    if (kjidim > 9000000) {
	s_stop("SOLVA_ERROR: max cubes exceeded", (ftnlen)31);
    }

/*     -- Prepare upto ncube cubes each containing upto nac atoms. The cube index */
/*     -- is kji. The atom index for each cube is in itab */

    for (l = 1; l <= 9000000; ++l) {
	itab[l - 1] = 0;
    }
    i__1 = *nats;
    for (l = 1; l <= i__1; ++l) {
	i__ = (xyz[l + 50000] - xmin) / rmax + 1.f;
	j = (xyz[l + 100000] - ymin) / rmax;
	k = (xyz[l + 150000] - zmin) / rmax;
	kji = k * jidim + j * idim + i__;
	n = itab[kji - 1] + 1;
	if (n > 150) {
	    s_stop("SOLVA_ERROR: max atoms per cube exceeded", (ftnlen)40);
	}
	itab[kji - 1] = n;
	natm[n + kji * 150 - 151] = l;
	cube[l - 1] = kji;
    }

/*     -- Process each atom in turn */

    nzp = 1.f / *zslice + .5f;
    i__1 = *nats;
    for (ir = 1; ir <= i__1; ++ir) {
	kji = cube[ir - 1];
	io = 0;
	area = 0.f;
	xr = xyz[ir + 50000];
	yr = xyz[ir + 100000];
	zr = xyz[ir + 150000];
	rr = rad[ir - 1];
	rrx2 = rr * 2.f;
	rrsq = radsq[ir - 1];

/*     -- Find the 'mkji' cubes neighboring the kji cube */

	for (k = -1; k <= 1; ++k) {
	    for (j = -1; j <= 1; ++j) {
		for (i__ = -1; i__ <= 1; ++i__) {
		    mkji = kji + k * jidim + j * idim + i__;
		    if (mkji >= 1) {
			if (mkji > kjidim) {
			    goto L14;
			}
			nm = itab[mkji - 1];
			if (nm >= 1) {

/*     -- record the atoms in inov that neighbor atom ir */

			    i__2 = nm;
			    for (m = 1; m <= i__2; ++m) {
				in = natm[m + mkji * 150 - 151];
				if (in != ir) {
				    ++io;
				    if (io > ict) {
					s_stop("SOLVA_ERROR: intrsctns > max",
						 (ftnlen)28);
				    }
				    dx[io - 1] = xr - xyz[in + 50000];
				    dy[io - 1] = yr - xyz[in + 100000];
/* Computing 2nd power */
				    r__1 = dx[io - 1];
/* Computing 2nd power */
				    r__2 = dy[io - 1];
				    dsq[io - 1] = r__1 * r__1 + r__2 * r__2;
				    d__[io - 1] = sqrt(dsq[io - 1]);
				    inov[io - 1] = in;
				}
			    }
			}
		    }
		}
	    }
	}
L14:
	if (io >= 1) {

/*     z resolution determined */

	    zres = rrx2 / nzp;
	    zgrid = xyz[ir + 150000] - rr - zres / 2.f;
	} else {
	    area = pix2 * rrx2;
	    goto L18;
	}

/*     section atom spheres perpendicular to the z axis */

	i__2 = nzp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    zgrid += zres;

/*     find the radius of the circle of intersection of */
/*     the ir sphere on the current z-plane */

/* Computing 2nd power */
	    r__1 = zgrid - zr;
	    rsec2r = rrsq - r__1 * r__1;
	    rsecr = sqrt(rsec2r);
	    i__3 = karc;
	    for (k = 1; k <= i__3; ++k) {
		arci[k - 1] = 0.f;
	    }
	    karc = 0;
	    i__3 = io;
	    for (j = 1; j <= i__3; ++j) {
		in = inov[j - 1];

/*     find radius of circle locus */

/* Computing 2nd power */
		r__1 = zgrid - xyz[in + 150000];
		rsec2n = radsq[in - 1] - r__1 * r__1;
		if (rsec2n <= 0.f) {
		    goto L10;
		}
		rsecn = sqrt(rsec2n);

/*     find intersections of n.circles with ir circles in section */

		if (d__[j - 1] >= rsecr + rsecn) {
		    goto L10;
		}

/*     do the circles intersect, or is one circle completely inside the other? */

		b = rsecr - rsecn;
		if (d__[j - 1] > dabs(b)) {
		    goto L20;
		}
		if (b <= 0.f) {
		    goto L9;
		}
		goto L10;

/*     if the circles intersect, find the points of intersection */

L20:
		++karc;
		if (karc >= ict) {
		    s_stop("SOLVA_ERROR: max intersections exceeded2", (
			    ftnlen)40);
		}

/*     Initial and final arc endpoints are found for the ir circle intersected */
/*     by a neighboring circle contained in the same plane. The initial endpoint */
/*     of the enclosed arc is stored in arci, and the final arc in arcf */
/*     law of cosines */

		trig_test__ = (dsq[j - 1] + rsec2r - rsec2n) / (d__[j - 1] * 
			2.f * rsecr);
		if (trig_test__ >= 1.f) {
		    trig_test__ = .99999f;
		}
		if (trig_test__ <= -1.f) {
		    trig_test__ = -.99999f;
		}
		alpha = acos(trig_test__);

/*     alpha is the angle between a line containing a point of intersection and */
/*     the reference circle center and the line containing both circle centers */

		beta = atan2(dy[j - 1], dx[j - 1]) + pi;

/*     beta is the angle between the line containing both circle centers and the x-axis */

		ti = beta - alpha;
		tf = beta + alpha;
		if (ti < 0.f) {
		    ti += pix2;
		}
		if (tf > pix2) {
		    tf -= pix2;
		}
		arci[karc - 1] = ti;
		if (tf >= ti) {
		    goto L3;
		}

/*     if the arc crosses zero, then it is broken into two segments. */
/*     the first ends at pix2 and the second begins at zero */

		arcf[karc - 1] = pix2;
		++karc;
L3:
		arcf[karc - 1] = tf;
L10:
		;
	    }

/*     find the accessible surface area for the sphere ir on this section */

	    if (karc != 0) {
		goto L19;
	    }
	    arcsum = pix2;
	    goto L25;

/*     The arc endpoints are sorted on the value of the initial arc endpoint */

L19:
	    sortag_(arci, &karc, tag);

/* *************************************** */
/*     calculate the accessible area */
/* *************************************** */

	    arcsum = arci[0];
	    t = arcf[tag[0] - 1];
	    if (karc == 1) {
		goto L11;
	    }
	    i__3 = karc;
	    for (k = 2; k <= i__3; ++k) {
		if (t < arci[k - 1]) {
		    arcsum = arcsum + arci[k - 1] - t;
		}
		tt = arcf[tag[k - 1] - 1];
		if (tt > t) {
		    t = tt;
		}
	    }
L11:
	    arcsum = arcsum + pix2 - t;

/*     The area/radius is equal to the accessible arc length x the section thickness. */

L25:
	    parea = arcsum * zres;

/*     Add the accessible area for this atom in this section to the area for this */
/*     atom for all the section encountered thus far */

	    area += parea;
L9:
	    ;
	}

/*     scale area to vdw shell */

L18:
	b = area * rr;
	accs[ir] = b;
/* ------------------------------------------------------------------ */
/* The following line converts from accessible to contact surface */
/*         c=(b*(rad(ir)-probe)**2)/(rad(ir)**2) */
/* ------------------------------------------------------------------ */
    }
    s_wsfe(&io___186);
    do_fio(&c__1, " SOLVA: PROGRAM ENDS CORRECTLY", (ftnlen)30);
    e_wsfe();
    return 0;
} /* solva_ */

/* Subroutine */ int sortag_(real *a, integer *n, integer *tag)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real t;
    static integer ij, il[16], tg, iu[16];
    static real tt;

    /* Parameter adjustments */
    --tag;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tag[i__] = i__;
    }
    m = 1;
    i__ = 1;
    j = *n;
L5:
    if (i__ >= j) {
	goto L70;
    }
L10:
    k = i__;
    ij = (j + i__) / 2;
    t = a[ij];
    if (a[i__] <= t) {
	goto L20;
    }
    a[ij] = a[i__];
    a[i__] = t;
    t = a[ij];
    tg = tag[ij];
    tag[ij] = tag[i__];
    tag[i__] = tg;
L20:
    l = j;
    if (a[j] >= t) {
	goto L40;
    }
    a[ij] = a[j];
    a[j] = t;
    t = a[ij];
    tg = tag[ij];
    tag[ij] = tag[j];
    tag[j] = tg;
    if (a[i__] <= t) {
	goto L40;
    }
    a[ij] = a[i__];
    a[i__] = t;
    t = a[ij];
    tg = tag[ij];
    tag[ij] = tag[i__];
    tag[i__] = tg;
    goto L40;
L30:
    a[l] = a[k];
    a[k] = tt;
    tg = tag[l];
    tag[l] = tag[k];
    tag[k] = tg;
L40:
    --l;
    if (a[l] > t) {
	goto L40;
    }
    tt = a[l];
L50:
    ++k;
    if (a[k] < t) {
	goto L50;
    }
    if (k <= l) {
	goto L30;
    }
    if (l - i__ <= j - k) {
	goto L60;
    }
    il[m - 1] = i__;
    iu[m - 1] = l;
    i__ = k;
    ++m;
    goto L80;
L60:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    ++m;
    goto L80;
L70:
    --m;
    if (m == 0) {
	return 0;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
L80:
    if (j - i__ >= 1) {
	goto L10;
    }
    if (i__ == 1) {
	goto L5;
    }
    --i__;
L90:
    ++i__;
    if (i__ == j) {
	goto L70;
    }
    t = a[i__ + 1];
    if (a[i__] <= t) {
	goto L90;
    }
    tg = tag[i__ + 1];
    k = i__;
L100:
    a[k + 1] = a[k];
    tag[k + 1] = tag[k];
    --k;
    if (t < a[k]) {
	goto L100;
    }
    a[k + 1] = t;
    tag[k + 1] = tg;
    goto L90;
} /* sortag_ */

/* Subroutine */ int summer_(char *sname, integer *slen, integer *nats, real *
	accs, integer *backbone, integer *polstats, integer *rtype, integer *
	resindex, char *resnam, real *ressums, real *tsums, ftnlen sname_len, 
	ftnlen resnam_len)
{
    /* System generated locals */
    integer i__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_rsfi(icilist *), e_rsfi(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, j;
    static real standarea[700]	/* was [100][7] */;
    static logical ok;
    static integer ir;
    extern integer readstring_(integer *, char *, integer *, ftnlen);
    static char res[3], line[256];
    static integer ilen, ires;
    static char acids[3*100];
    extern integer fopen_(integer *, char *, integer *, char *, ftnlen, 
	    ftnlen);
    static logical stand;
    extern /* Subroutine */ int which3_(char *, char *, integer *, integer *, 
	    logical *, ftnlen, ftnlen);
    static integer nacids, rindex[5000];

    /* Fortran I/O blocks */
    static cilist io___199 = { 0, 3, 0, "(4a)", 0 };
    static icilist io___204 = { 0, line+16, 0, "(f7.2)", 7, 1 };
    static icilist io___206 = { 0, line+29, 0, "(f7.2)", 7, 1 };
    static icilist io___207 = { 0, line+42, 0, "(f7.2)", 7, 1 };
    static icilist io___208 = { 0, line+55, 0, "(f7.2)", 7, 1 };
    static icilist io___209 = { 0, line+68, 0, "(f7.2)", 7, 1 };
    static icilist io___210 = { 0, line+81, 0, "(f7.2)", 7, 1 };
    static icilist io___211 = { 0, line+94, 0, "(f7.2)", 14, 1 };
    static cilist io___212 = { 0, 4, 0, "(a,i3,a)", 0 };
    static cilist io___213 = { 0, 4, 0, "(a)", 0 };


/*     -- 	program to sum atomic accessibilities by residue. */
/*     --	copes with atom and hetatom records */
/*     --	produces relative accessibilities for the 20 common aminos */
/*     --	ouput written to .rsa file (channel 4) */

/* Parameters uses in ACCALL.F */
/*     --  PARAMETERS: maxa = max atoms per residue */
/*     --	      maxr = max number of "standard" residue types */
/*     --              maxs = max number of atoms in PDB file */
/*     --              maxx = max number of residues in PDB file */
/*     --              maxc = max number of chains */

/* ============================================================= */

/* SOLVA LOCAL PARAMETERS: */
/* ncube = maximum number of cubes allowed for placing of atoms */
/* nac   = maximum number of atoms per cube */
/* nint  = maximum number of sphere intersections */


/*     if "standard.data" exists in, read them in. */

    /* Parameter adjustments */
    --tsums;
    ressums -= 40001;
    resnam -= 10;
    --resindex;
    --rtype;
    --polstats;
    --backbone;
    --accs;

    /* Function Body */
    stand = FALSE_;
    if (fopen_(&c__1, sname, slen, "old", (ftnlen)256, (ftnlen)3) != 0) {
	s_wsfe(&io___199);
	do_fio(&c__1, "REM  Relative accessibilites read from", (ftnlen)38);
	do_fio(&c__1, " external file \"", (ftnlen)16);
	do_fio(&c__1, sname, (*slen));
	do_fio(&c__1, "\"", (ftnlen)1);
	e_wsfe();
	stand = TRUE_;
	i__ = 0;
	while(readstring_(&c__1, line, &ilen, (ftnlen)256) >= 0 && i__ < 100) 
		{
	    if (s_cmp(line, "ATOM", (ftnlen)4, (ftnlen)4) == 0) {
		++i__;
		s_copy(acids + (i__ - 1) * 3, line + 12, (ftnlen)3, (ftnlen)3)
			;
		s_rsfi(&io___204);
		do_fio(&c__1, (char *)&standarea[i__ - 1], (ftnlen)sizeof(
			real));
		e_rsfi();
		s_rsfi(&io___206);
		do_fio(&c__1, (char *)&standarea[i__ + 99], (ftnlen)sizeof(
			real));
		e_rsfi();
		s_rsfi(&io___207);
		do_fio(&c__1, (char *)&standarea[i__ + 199], (ftnlen)sizeof(
			real));
		e_rsfi();
		s_rsfi(&io___208);
		do_fio(&c__1, (char *)&standarea[i__ + 299], (ftnlen)sizeof(
			real));
		e_rsfi();
		s_rsfi(&io___209);
		do_fio(&c__1, (char *)&standarea[i__ + 399], (ftnlen)sizeof(
			real));
		e_rsfi();
		s_rsfi(&io___210);
		do_fio(&c__1, (char *)&standarea[i__ + 499], (ftnlen)sizeof(
			real));
		e_rsfi();
		s_rsfi(&io___211);
		do_fio(&c__1, (char *)&standarea[i__ + 599], (ftnlen)sizeof(
			real));
		e_rsfi();
	    }
	}
	s_wsfe(&io___212);
	do_fio(&c__1, " RELATIVE (STANDARD) ACCESSIBILITIES READFOR ", (
		ftnlen)45);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, " AMINO ACIDS", (ftnlen)12);
	e_wsfe();
	cl__1.cerr = 0;
	cl__1.cunit = 1;
	cl__1.csta = 0;
	f_clos(&cl__1);
    } else {
	s_wsfe(&io___213);
	do_fio(&c__1, " NO STANDARD VALUES INPUT", (ftnlen)25);
	e_wsfe();
    }
    nacids = i__;
    i__1 = resindex[*nats];
    for (i__ = 1; i__ <= i__1; ++i__) {
	rindex[i__ - 1] = 0;
	if (stand) {
	    s_copy(res, resnam + i__ * 10, (ftnlen)3, (ftnlen)3);
	    which3_(res, acids, &ires, &nacids, &ok, (ftnlen)3, (ftnlen)3);
	    if (ok) {
		rindex[i__ - 1] = ires;
	    }
	}
    }

/*     -- sum the values */

    i__1 = *nats;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ir = resindex[i__];
	tsums[1] += accs[i__];
	ressums[ir + 40000] += accs[i__];
	if (backbone[i__] == 0) {
	    ressums[ir + 60000] += accs[i__];
	    tsums[5] += accs[i__];
	} else {
	    ressums[ir + 55000] += accs[i__];
	    tsums[4] += accs[i__];
	    if (polstats[i__] == 0) {
		ressums[ir + 45000] += accs[i__];
		tsums[2] += accs[i__];
	    } else if (polstats[i__] == 1) {
		ressums[ir + 50000] += accs[i__];
		tsums[3] += accs[i__];
	    }
	}
	if (polstats[i__] == 0) {
	    ressums[ir + 65000] += accs[i__];
	    tsums[6] += accs[i__];
	} else {
	    ressums[ir + 70000] += accs[i__];
	    tsums[7] += accs[i__];
	}
    }

/*     -- calculate relative accessibilities */

    i__1 = resindex[*nats];
    for (i__ = 1; i__ <= i__1; ++i__) {
	ires = rindex[i__ - 1];
	if (stand && ires != 0) {
	    for (j = 1; j <= 7; ++j) {
		if (standarea[ires + j * 100 - 101] > 0.f) {
		    ressums[i__ + (j + 14) * 5000] = ressums[i__ + (j + 7) * 
			    5000] * 100.f / standarea[ires + j * 100 - 101];
		} else {
		    ressums[i__ + (j + 14) * 5000] = 0.f;
		}
	    }
	} else {
	    for (j = 1; j <= 7; ++j) {
		ressums[i__ + (j + 14) * 5000] = -99.9f;
	    }
	}
    }
    return 0;
} /* summer_ */

/* Subroutine */ int polguess_(char *atom, integer *ip, ftnlen atom_len)
{
    *ip = 0;
    if (*(unsigned char *)&atom[1] == 'O' || *(unsigned char *)&atom[1] == 
	    'N' || *(unsigned char *)&atom[1] == 'A') {
	*ip = 1;
    }
    return 0;
} /* polguess_ */

integer what_atom__(char *atom, integer *ir, integer *n, ftnlen atom_len)
{
    /* Initialized data */

    static char mc[4*5] = " N  " " C  " " O  " " OXT" " CA ";
    static char nc[4*11] = " P  " " O1P" " O2P" " O5*" " C5*" " C4*" " O4*" 
	    " C3*" " O3*" " C2*" " C1*";

    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;

    ret_val = 0;
    if (*ir == 1) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (s_cmp(atom, mc + (i__ - 1 << 2), (ftnlen)4, (ftnlen)4) == 0) {
		return ret_val;
	    }
	}
    } else if (*ir == 2) {
	for (i__ = 1; i__ <= 11; ++i__) {
	    if (s_cmp(atom, nc + (i__ - 1 << 2), (ftnlen)4, (ftnlen)4) == 0) {
		return ret_val;
	    }
	}
    } else {
	for (i__ = 1; i__ <= 4; ++i__) {
	    if (s_cmp(atom, mc + (i__ - 1 << 2), (ftnlen)4, (ftnlen)4) == 0) {
		return ret_val;
	    }
	}
	for (i__ = 1; i__ <= 11; ++i__) {
	    if (s_cmp(atom, nc + (i__ - 1 << 2), (ftnlen)4, (ftnlen)4) == 0) {
		return ret_val;
	    }
	}
    }
    ret_val = 1;
    return ret_val;
} /* what_atom__ */

/* Subroutine */ int vanin_(char *vname, integer *vlen, integer *nacids, char 
	*aacids, char *anames, integer *numats, real *vradii, integer *spolar,
	 integer *rtype, ftnlen vname_len, ftnlen aacids_len, ftnlen 
	anames_len)
{
    /* System generated locals */
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_clos(cllist *);

    /* Local variables */
    static char c__[256*256];
    static integer i__, l[256], n;
    extern doublereal readfloat_(char *, integer *, ftnlen);
    static char aa3[3];
    extern integer readstring_(integer *, char *, integer *, ftnlen);
    static char card[256];
    static integer ilen, nlabs;
    extern integer fopen_(integer *, char *, integer *, char *, ftnlen, 
	    ftnlen), parse_(char *, integer *, char *, char *, integer *, 
	    ftnlen, ftnlen, ftnlen), readint_(char *, integer *, ftnlen);


/* -- Read in van der Waal radii from external file "vfile" */
/* -- nacids = number of residues read in */
/* -- aacids = *3 character array containing amino acid names */
/* -- anames = *4 character array containing atom names for each residue */
/* -- numats = array containing number of atoms for each residue */
/* -- radii  = vdw radii for each atom, indexed identically to anames */


/* Parameters uses in ACCALL.F */
/*     --  PARAMETERS: maxa = max atoms per residue */
/*     --	      maxr = max number of "standard" residue types */
/*     --              maxs = max number of atoms in PDB file */
/*     --              maxx = max number of residues in PDB file */
/*     --              maxc = max number of chains */

/* ============================================================= */

/* SOLVA LOCAL PARAMETERS: */
/* ncube = maximum number of cubes allowed for placing of atoms */
/* nac   = maximum number of atoms per cube */
/* nint  = maximum number of sphere intersections */

    /* Parameter adjustments */
    --rtype;
    spolar -= 101;
    vradii -= 101;
    --numats;
    anames -= 404;
    aacids -= 3;

    /* Function Body */
    if (fopen_(&c__10, vname, vlen, "old", (ftnlen)256, (ftnlen)3) == 0) {
	s_stop("ERROR: unable to open \"vdw.radii\"", (ftnlen)33);
    }
    *nacids = 0;
    nlabs = 0;
    while(readstring_(&c__10, card, &ilen, (ftnlen)256) >= 0) {
	n = parse_(card, &ilen, " ", c__, l, (ftnlen)256, (ftnlen)1, (ftnlen)
		256);
	if (s_cmp(c__, "RESIDUE", (ftnlen)256, (ftnlen)7) == 0) {
	    ++(*nacids);
	    rtype[*nacids] = 1;
	    if (s_cmp(c__ + 256, "NUCL", (ftnlen)4, (ftnlen)4) == 0) {
		rtype[*nacids] = 2;
	    }
	    if (s_cmp(c__ + 256, "HETA", (ftnlen)4, (ftnlen)4) == 0) {
		rtype[*nacids] = 3;
	    }
	    if (*nacids > 100) {
		s_stop("ERROR: increase maxr", (ftnlen)20);
	    }
	    s_copy(aa3, c__ + 512, (ftnlen)3, (ftnlen)3);
	    for (i__ = 1; i__ <= 3; ++i__) {
		if (*(unsigned char *)&aa3[i__ - 1] == '_') {
		    *(unsigned char *)&aa3[i__ - 1] = ' ';
		}
	    }
	    s_copy(aacids + *nacids * 3, aa3, (ftnlen)3, (ftnlen)3);
	    numats[*nacids] = 0;
	}
	if (s_cmp(c__, "ATOM", (ftnlen)256, (ftnlen)4) == 0) {
	    ++numats[*nacids];
	    if (numats[*nacids] > 100) {
		s_stop("ERROR: increase maxa", (ftnlen)20);
	    }
/*            aa4 = card(6:9) */
/*            do i = 1, 4 */
/*               if(aa4(i:i).eq.'_')aa4(i:i)=' ' */
/*            enddo */
	    s_copy(anames + (*nacids + numats[*nacids] * 100 << 2), card + 5, 
		    (ftnlen)4, (ftnlen)4);
	    vradii[*nacids + numats[*nacids] * 100] = readfloat_(card + 10, &
		    c__4, (ftnlen)4);
/*            write(6,'(2a,1x,a,f8.3)') */
/*     -           'AT> ',anames(nacids,numats(nacids)), */
/*     -           aacids(nacids),vradii(nacids,numats(nacids)) */
	    if (n >= 4) {
		spolar[*nacids + numats[*nacids] * 100] = readint_(card + 15, 
			&c__1, (ftnlen)1);
	    } else {
		spolar[*nacids + numats[*nacids] * 100] = -1;
	    }
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* vanin_ */

/* Main program alias */ int accall_ () { MAIN__ (); return 0; }
