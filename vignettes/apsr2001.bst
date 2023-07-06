% BibTeX bibliography style `apsr2001' 
% Version 1.01
% 2006 September 19
%
% ACKNOWLEDGEMENT
%   This document is a modified version of Peter Williams' apsr.bst,
%   which was in turn based on the original alpha.bst.
%
% MODIFIED BY
%   John G. Bullock, john.bullock@aya.yale.edu 
%
% INSTALLATION INSTRUCTIONS
%   Find your local copy of apsr.bst.  (On my hard drive, it's in
%   c:/texmf/bibtex/bst/harvard/.)  Put this file in the same
%   directory.  Update your file name database.  (If you use MiKTeX,
%   this just entails opening "MiKTeX Options" and clicking "Refresh 
%   Now".)
%        In the preambles to your LaTeX files, specify 
%
%        \bibliographystyle{apsr2001}
%
%   instead of 
%
%        \bibliographystyle{apsr2001}
%
% REASON FOR EXISTING
%   Current APSR style rules are set out in the 2001 edition 
%   of the APSA's Style Manual for Political Science.  The 
%   original apsr.bst goes a long way toward meeting those 
%   rules.  But there are at least 18 ways in which it falls 
%   short and produces rather ugly bibliographies.  This file
%   fixes most of the problems.  See the changelog for details.
%
% CHANGELOG
%   Version 1.01, 2006 September 19:
%   --Restored en dashes between page numbers of cited articles.
%
%   Version 1.0,  2006 May 5:
%   --Added 'translator' field for books.  I don't know how robust
%     it is, i.e., how it will work with complicated bibliography
%     entries of the "book" type.  For multiple translators, I may
%     want to borrow the code in format.authors and format.rem.names.
%
%   --When an article cites two articles written by the same 
%     author in the same year, a letter is typically appended 
%     to distinguish them, e.g., Converse 1964a vs. Converse 1964b.
%     apsr.bst italicized these appended letters.  apsr2001.bst 
%     does not: I edited {reverse.pass} to prevent the needless 
%     italicization.
%
%   --Changed format.articles and spacing in format.vol.num.pages 
%     so wire reports from news agencies (e.g., the Associated Press) 
%     will format OK.  
%
%   --Often, we want to cite chapters that appear in edited volumes.
%     Typically, the chapters have the INCOLLECTION or INPROCEEDINGS
%     entry type.  If more than one chapter is to be cited from the 
%     same volume, the most efficient thing to do is to cross-reference
%     (with the CROSSREF directive) a PROCEEDINGS entry that contains 
%     all of the information common to all chapters: the book title, 
%     editor(s) of the collection, publication year, and so forth.  But 
%     cross-referencing wasn't working.  Now it is.
%
%   --Changed the proceedings function, which serves for edited volumes.
%     It now takes "booktitle" as well as "title."  If either is missing,
%     a BibTeX warning will be generated; like many BibTeX warnings, it
%     can be ignored.  This was done to make cross-referencing work; see
%     above.
%
%   --Changed the conferences function, which is also the inproceedings
%     function, to make it more workable for real conference papers.
%     (See APSA Style Manual 2001, 32.)
%
%   --Per the Style Manual, stopped n.dashifying the dashes between 
%     page numbers.  
%
%   --Changed format.vol.num.pages to degrade gracefully, as required by
%     APSA Style Manual 2001, 31.  Uses journal month/season if available, 
%     followed by journal number.  Also added a space between the 
%     volume number and the opening parenthesis of the month/season/issue
%     number; this, too, is in keeping with the Manual.
%
%   --If a journal number is printed in the bibliography, the space
%     between the volume number and the journal number will be
%     non-breaking, so that no line in the bibliography will begin 
%     with a journal number.  E.g., an entry will never look like
%
%          Smith, John.  "The Article."  American Political Science Review 10 
%               (2): 147-57.
%
%   --Added format.rem.names.editors, which is based on format.rem.names.  
%     Also changed format.editors and format.editors.reverse to refer
%     to the new function rather than to format.rem.names.
%         This was done so that, in bibliographies, incollection
%     entries with two editors would appear with no comma (e.g., 
%     "ed. John Smith and Jane Doe") while entries with three or more
%     editors would have the appropriate final comma (e.g., 
%     "ed. John Smith, Jane Doe, and Paul Sniderman").
%
%   --Changed format.vol.num.pages to allow for a space between the 
%     colon and the page numbers.
%
%   --Changed FUNCTION {incollection} so that chapter pages aren't 
%     entered into bibliography.  (APSA style guide 2001, 29-30: 
%         pages are not required.)
%   
%   --Changed FUNCTION {phdthesis} to accord with APSA style guide:
%     Removed location printing and changed "PhD thesis" to "Ph.D. 
%     diss."
%
%   --In format.lab.names.full, changed
%               { " \harvardand\ " * t * } 
%			 to
%	    	      { numnames #2 >	
%	        		{ ", \harvardand\ " * t * }
%	  	          { " \harvardand\ " * t * }
%			   if$ 
% 			 }	
%
%   --Added ' ' after 'In' at format.in.ed.booktitle.
%     In format.rem.names, added ", " to the "harvard and" line.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    originalyear
    pages
    publisher
    school
    series
    title
    translator
    type
    URL
    volume
    year
  }
  { field.used etal.allowed etal.required} %%%XXX change
  { extra.label sort.label list.year }

INTEGERS { output.state before.all mid.sentence after.sentence after.block }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}

STRINGS { s t f originalyearstring}

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { " " * write$ }
    { output.state after.block =
    { add.period$ write$
      newline$
      "\newblock " write$
    }
    { output.state before.all =
        'write$
        { add.period$ " " * write$ }
      if$
    }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {item.check}
{ 't :=
  empty$
    { "empty " t * " in " * cite$ * warning$ }
    { skip$ }
  if$
}

FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{\em " swap$ * "}" * }
  if$
}

FUNCTION {embolden}
{ duplicate$ empty$
    { pop$ "" }
    { "{\bf " swap$ * "}" * }
  if$
}

FUNCTION {quote}
{ duplicate$ empty$
    { pop$ "" }
    { "``" swap$ * "''" * }
  if$
}

FUNCTION {write.url}
{ URL empty$
    { skip$ }
    { "\newline\harvardurl{" URL * "}" * write$ newline$ }
  if$
}

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.rem.names.editors} {
  's :=
  'f :=
  #2 'nameptr :=
  s num.names$ 'numnames :=
  numnames #1 - 'namesleft :=
  { namesleft #0 > }
    { s nameptr f format.name$ 't :=
      nameptr #1 >
        { namesleft #1 >
            { ", " * t * }
            { t "others" =
                { " et~al." * }
                { % begin if loop for editor commas in citations
                  numnames #2 >
                    { ", \harvardand\ " * t * }
                    { " \harvardand\ " * t * }
                  if$
                } % end if loop for editor commas in citations
            if$
            }
        if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.rem.names} {
   's :=
   'f :=
   #2 'nameptr :=
   s num.names$ 'numnames :=
   numnames #1 - 'namesleft :=
   { namesleft #0 > }
     { s nameptr f format.name$ 't :=
       nameptr #1 >
         { namesleft #1 >
             { ", " * t * }
             { t "others" =
                 { " et~al." * }
                 { ", \harvardand\ " * t * }
             if$
             }
         if$
         }
         't
       if$
       nameptr #1 + 'nameptr :=
       namesleft #1 - 'namesleft :=
    }
   while$
 }


FUNCTION {format.authors}
{ author empty$
    { "" }
    { author #1 "{vv~}{ll}{, jj}{, ff}" format.name$
      "{ff }{vv~}{ll}{, jj}" author format.rem.names
    }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor #1 "{vv~}{ll}{, jj}{, ff}" format.name$
      "{ff }{vv~}{ll}{, jj}" editor format.rem.names.editors
      editor num.names$ #1 >
    { ", eds" * }
    { ", ed." * }
      if$
    }
  if$
}

FUNCTION {format.editors.reverse}
{ editor empty$
    { "" }
    { ", ed. "
      editor #1 "{ff }{vv~}{ll}{, jj}" format.name$ *
      "{ff }{vv~}{ll}{, jj}" editor format.rem.names.editors
    }
  if$
}

FUNCTION {format.originalyear}
{ originalyear empty$
	{ "" }                      		% if originalyear field is empty
	{ "[" originalyear * "] " * }      % if it isn't empty
	if$
}

FUNCTION {format.title}
{ title empty$
    { "" }
    { title add.period$ }
  if$
}

FUNCTION {format.translator}
{ translator empty$
	{ "" }                      		% if translator field is empty
	{ "Trans. " translator * ". " * }   % if it isn't empty
	if$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
    { t #1 #2 substring$ "--" = not
        { "--" *
          t #2 global.max$ substring$ 't :=
        }
        {   { t #1 #1 substring$ "-" = }
        { "-" *
          t #2 global.max$ substring$ 't :=
        }
          while$
        }
      if$
    }
    { t #1 #1 substring$ *
      t #2 global.max$ substring$ 't :=
    }
      if$
    }
  while$
}

FUNCTION {format.btitle}
{ title emphasize
}

FUNCTION {format.booktitle}
{ booktitle emphasize
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { "Vol." volume tie.or.space.connect
      series empty$
    'skip$
    { " of " * series emphasize * }
      if$
      "volume and number" number either.or.check
    }
  if$
}

FUNCTION {format.number.series}
{ volume empty$
    { number empty$
    { series field.or.null }
    { output.state mid.sentence =
        { "number" }
        { "Number" }
      if$
      number tie.or.space.connect
      series empty$
        { "there's a number but no series in " cite$ * warning$ }
        { " {\em in} " * series quote * }
      if$
    }
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
    { edition "l" change.case$ " ed." * } 
    { edition "t" change.case$ " ed." * } 
      if$
    }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
    { #1 'multiresult := }
    { t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
    { "pp.~" pages n.dashify * }
    { "pp.~" pages * }
    { "p.~" pages * }
      if$
    }
  if$
}

FUNCTION {format.address}
{ address empty$
    { "" }
    { address ": " *
    }
  if$
}

FUNCTION {format.publisher}
{ publisher empty$
    { "" }
    { publisher }    
  if$
}

FUNCTION {format.vol.num.pages}
%{ volume field.or.null
{ volume empty$
	{""}
	{ " " * volume field.or.null }
  if$
  month empty$
    { number empty$
        'skip$ % do nothing if "number" field is empty
        { "~(" number * ")" * *                             % { " (" number * ")" * *
        volume empty$
            { "there's a number but no volume in " cite$ * warning$ }
            'skip$
        if$
        }
    if$ }
    { " (" month * ")" * * }
  if$
  pages empty$
    'skip$
    { duplicate$ empty$
      { pop$ " " * format.pages }      % { pop$ format.pages }      
	 { ": " * pages * }
    if$
    }
  if$
}


FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
    { "chapter" }
    { type "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
    'skip$
    { ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
    { " In " booktitle emphasize * }
    { " In " booktitle emphasize * format.editors.reverse * }
      if$
    }
  if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
  key empty$ not and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Technical Report" }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}

FUNCTION {format.article.crossref}
{ key empty$
    { journal empty$
    { "need key or journal for " cite$ * " to crossref " * crossref *
      warning$
      ""
    }
    { "in {\em " journal * "\/} \cite{" * crossref * "}" *}
      if$
    }
    { "{\em in} \citeasnoun{" crossref * "}" * }
  if$
 
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "in "
    }
    { "Vol." volume tie.or.space.connect
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
    { series empty$
        { "need editor, key, or series for " cite$ * " to crossref " *
          crossref * warning$
          "" *
        }
        { "{\em " * series * "\/} \cite{" * crossref * "}" *}
      if$
    }
    { " \citeasnoun{" * crossref * "}" * }
      if$
    }
    { " \citeasnoun{" * crossref * "}" * }
  if$
}

FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
    { booktitle empty$
        { "need editor, key, or booktitle for " cite$ * " to crossref " *
          crossref * warning$
          ""
        }
        { "in {\em " booktitle * "\/}" * " \cite{" * crossref * "}" *}
      if$
    }
    { "{\em in} \citeasnoun{" crossref * "}" * }
      if$
    }
    %{ "In {\em " booktitle * ",\/} ed. " format.editors * " \citeasnoun{" crossref * "}" * }
    { " In " * format.booktitle *
	 format.editors.reverse * ". " * 
	 format.address *
	 format.publisher 
      % publisher "publisher" output.check 
    }
  if$
  
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}

INTEGERS { ind tsslen }

STRINGS { tss ret rss istr }

FUNCTION {replace.substring}{
  'rss :=
  'tss :=
  'istr :=
  "" 'ret :=
  tss text.length$ 'tsslen :=
  #1 'ind :=
    { istr ind tsslen substring$ "" = not }
    { istr ind tsslen substring$ tss =
        { ret rss * 'ret :=
          ind tsslen + 'ind :=
        }
        { ret istr ind #1 substring$ * 'ret :=
          ind #1 + 'ind :=
        }
      if$
    }
  while$
  ret
}

FUNCTION {format.lab.names.full} 
{
 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr "{vv~}{ll}" format.name$ 't :=
      nameptr #1 >
    { namesleft #1 >
        { ", " * t * }
        { t "others" =
        { " et~al." * }
	   { numnames #2 >	
	        { ", \harvardand\ " * t * }
	        { " \harvardand\ " * t * }
		if$ 
	   }	
          if$
        }
      if$
    }
    't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}




FUNCTION {format.lab.names.abbr}
{ 's :=
  s num.names$ 'numnames :=
  numnames #1 >
    { numnames #3 >
    { s #1 "{vv~}{ll}" format.name$ " et~al." * }
    { s format.lab.names.full }
      if$
    }
    { s #1 "{vv~}{ll}" format.name$ }
  if$
}

INTEGERS { author.field editor.field organization.field title.field key.field booktitle.field }

FUNCTION {init.field.constants}
{ #0 'author.field :=
  #1 'editor.field :=
  #2 'organization.field :=
  #3 'title.field :=
  #4 'key.field :=
  #5 'booktitle.field :=
}

FUNCTION {make.list.label}
{ author.field field.used =
    { format.authors }
    { editor.field field.used =
        { format.editors }
        { organization.field field.used =
            { "The " #4 organization chop.word #3 text.prefix$ }
            { title.field field.used =
                { format.btitle }
                { key.field field.used =
                    { key #3 text.prefix$ }
                    { "Internal error :001 on " cite$ * " label" * warning$ }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
  if$
}

FUNCTION {make.full.label}
{ author.field field.used =
    { author format.lab.names.full }
    { editor.field field.used =
        { editor format.lab.names.full }
        { organization.field field.used =
            { "The " #4 organization chop.word #3 text.prefix$ }
            { title.field field.used =
                { format.btitle }
                { key.field field.used =
                    { key #3 text.prefix$ }
                    { "Internal error :001 on " cite$ * " label" * warning$ }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
  if$
}

FUNCTION {make.abbr.label} %%%XXX change
{
  etal.allowed
    { author.field field.used =
        { author format.lab.names.abbr }
        { editor.field field.used =
            { editor format.lab.names.abbr }
            { organization.field field.used =
                { "The " #4 organization chop.word #3 text.prefix$ }
                { title.field field.used =
                    { format.btitle }
                    { key.field field.used =
                        { key #3 text.prefix$ }
                        {"Internal error :001 on " cite$ * " label" * warning$ }
                      if$
                    }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
    { make.full.label }
  if$
}

FUNCTION {output.bibitem}
{ newline$
  etal.allowed        %%%XXX change
  etal.required
  and
    {
      "\harvarditem[" write$
      make.abbr.label write$
      "]{" write$
    }
    {
      "\harvarditem{" write$
    }
  if$
  make.full.label write$
  "}{" write$
  list.year write$
  "}{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {list.label.output}
{ make.list.label add.period$ " " * write$
}

FUNCTION {article}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  author "author" item.check
  title.field field.used =
    { skip$ }
    { format.title quote "title" output.check }
  if$
  crossref missing$
    { journal emphasize "journal" duplicate$ item.check
      year empty$
        { ". Forthcoming" * output new.block }
        { "" * format.vol.num.pages * output }  % { " " * format.vol.num.pages * output }
      if$
     }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  new.block
  note output
  fin.entry
  write.url
}

FUNCTION {book}
{ output.bibitem
  list.label.output
  %format.originalyear output
  list.year output.nonnull
  new.block
  author empty$
    { editor "author and editor" item.check }
    { crossref missing$
    { "author and editor" editor either.or.check }
    'skip$
      if$
    }
  if$
  title.field field.used =
    { skip$ }
    { format.btitle "title" output.check }
  if$
  new.block
  crossref missing$
    { format.bvolume output
	 format.translator output
      format.number.series output
      format.edition output
      format.address output
      publisher "publisher" output.check
    }
    { format.book.crossref output.nonnull
      format.edition output
    }
  if$
  new.block
  note output
  fin.entry
  write.url  
  %format.translator write$
  %newline$
}

FUNCTION {booklet}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.title quote "title" output.check }
  if$
  format.address output
  howpublished output
  new.block
  note output
  fin.entry
  write.url
}

FUNCTION {inbook}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  author empty$
    { editor "author and editor" item.check }
    { crossref missing$
    { "author and editor" editor either.or.check }
    'skip$
      if$
    }
  if$
  title.field field.used =
    { skip$ }
    { format.btitle "title" output.check }
  if$
  new.block
  crossref missing$
    { format.bvolume output
      format.number.series output
      format.edition output
      format.address output
      publisher "publisher" output.check
    }
    { format.book.crossref output.nonnull
      format.edition output
    }
  if$
  format.chapter.pages "chapter and pages" output.check
  new.block
  note output
  fin.entry
  write.url
}

FUNCTION {incollection}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.title quote "title" output.check }
  if$
  author "author" item.check
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      new.block
      format.edition output
      format.bvolume output
      format.number.series output
      format.address output
      publisher "publisher" output.check
    }
    { format.incoll.inproc.crossref output.nonnull
    }
  if$
  %format.chapter.pages output
  new.block
  note output
  fin.entry
  write.url
}


FUNCTION {inproceedings}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.title quote "title" output.check }
  if$
  author "author" item.check
  crossref missing$
    { booktitle empty$
        { skip$ }
	   { 
          format.in.ed.booktitle "booktitle" output.check
          new.block
        }
      if$
      format.bvolume output
      format.number.series output
      address empty$
        { organization output
          publisher output
        }
        { organization output
              %format.address output.nonnull
    	    	    ", " * address output.nonnull
          %publisher output
        }
      if$
    }
    { format.incoll.inproc.crossref output.nonnull
    }
  if$
  format.pages output
  new.block
  note output
  fin.entry
  write.url
} 

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.btitle "title" output.check }
  if$
  new.block
  format.edition output
  new.block
  author empty$
    { organization empty$
    { format.address output }
    'skip$
      if$
    }
    { format.address output
      organization output  
    }
  if$
  new.block
  note output
  fin.entry
  write.url
}

FUNCTION {mastersthesis}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  author "author" item.check
  title.field field.used =
    { skip$ }
    { format.title "title" output.check }
  if$
  "Master's thesis" format.thesis.type output.nonnull
  school "school" output.check
  format.address output
  new.block
  note output
  fin.entry
  write.url
}

FUNCTION {misc}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  title.field field.used =
    { skip$ }
    { format.title quote output }
  if$
  howpublished output
  new.block
  note output
  fin.entry
  write.url
  empty.misc.check
}


FUNCTION {phdthesis}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  author "author" item.check
  title.field field.used =
    { skip$ }
    { format.title quote "title" output.check }
  if$
  "Ph.D. diss." format.thesis.type output.nonnull
  school "school" output.check
  %format.address output
  new.block
  note output
  fin.entry
  write.url
}

FUNCTION {proceedings}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  %title.field field.used =
  %  { skip$ }
  %  { format.btitle "title" output.check }
  %if$
  booktitle.field field.used =
    { skip$ }
    { format.booktitle "booktitle" output.check }
  if$
  new.block
  format.bvolume output
  format.number.series output
  address empty$
    { editor empty$
    { skip$ }
    { organization output
    }
      if$
      publisher output
    }
    { editor empty$
    'skip$
    { organization output }
      if$
      format.address output.nonnull
      publisher output
    }
  if$
  new.block
  note output
  fin.entry
  write.url
}


FUNCTION {techreport}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  author "author" item.check
  title.field field.used =
    { skip$ }
    { format.title "title" output.check }
  if$
  format.tr.number output.nonnull
  institution "institution" output.check
  format.address output
  new.block
  note output
  fin.entry
  write.url
}

FUNCTION {unpublished}
{ output.bibitem
  list.label.output
  list.year output.nonnull
  new.block
  author "author" item.check
  title.field field.used =
    { skip$ }
    { format.title quote "title" output.check }
  if$
  note "note" output.check
  fin.entry
  write.url
}

FUNCTION {default.type} { misc }

MACRO {jan} {"January"}
MACRO {feb} {"February"}
MACRO {mar} {"March"}
MACRO {apr} {"April"}
MACRO {may} {"May"}
MACRO {jun} {"June"}
MACRO {jul} {"July"}
MACRO {aug} {"August"}
MACRO {sep} {"September"}
MACRO {oct} {"October"}
MACRO {nov} {"November"}
MACRO {dec} {"December"}
MACRO {tcs} {"Theoretical Computer Science"}

MACRO {APSR} {"American Political Science Review"}

READ

EXECUTE {init.field.constants}

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

FUNCTION {sortify.names}
{ " \harvardand\ " " " replace.substring
  " et~al." " zzz" replace.substring
  sortify
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
    { title.field 'field.used := }
    { key.field 'field.used := }
      if$
    }
    { author.field 'field.used := }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
    { key empty$
        { title.field 'field.used := }
        { key.field 'field.used := }
      if$
    }
    { editor.field 'field.used := }
      if$
    }
    { author.field 'field.used := }
  if$
}

FUNCTION {author.key.organization.label}
{ author empty$
    { key empty$
    { organization empty$
        { title.field 'field.used := }
        { organization.field 'field.used := }
      if$
    }
    { key.field 'field.used := }
      if$
    }
    { author.field 'field.used := }
  if$
}

FUNCTION {editor.key.organization.label}
{ editor empty$
    { key empty$
    { organization empty$
        { title.field 'field.used := }
        { organization.field 'field.used := }
      if$
    }
    { key.field 'field.used := }
      if$
    }
    { editor.field 'field.used := }
  if$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

% in btxbst.doc, determines what sort of in-text citation label 
% is made based on the entry type (e.g., book, article)
FUNCTION {calc.label} %%%XXX change 
{ 
  make.abbr.label 
  title.field field.used =
    { sort.format.title }
    { sortify.names }
  if$  
  year field.or.null purify$ #-1 #4 substring$ sortify
  *
  'sort.label :=
}

FUNCTION {preliminaries} %%%XXX change
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
    'editor.key.organization.label
    { type$ "manual" =
        'author.key.organization.label
        'author.key.label
      if$
    }
      if$
    }
  if$
  author.field field.used = %%%XXX change
    {
      author num.names$ #2 >
        { #1 }
        { #0 }
      if$
      'etal.required :=
    }
    {
      editor.field field.used = 
        {
          editor num.names$ #2 >
            { #1 }
            { #0 }
          if$
        }
        { #0 }
      if$
      'etal.required :=
    }
  if$
  #1 'etal.allowed :=
}

FUNCTION {first.presort}
{ calc.label
  sort.label
  title.field field.used =
    { skip$ }
    { "    "
      *
      make.list.label sortify.names
      *
      "    "
      *
      title field.or.null
      sort.format.title
      *
    }
  if$
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {preliminaries}

ITERATE {first.presort}

SORT

STRINGS { last.sort.label next.extra last.full.label}

INTEGERS { last.extra.num last.etal.allowed}

FUNCTION {initialize.confusion}
{ #0 int.to.chr$ 'last.sort.label :=
  #0 int.to.chr$ 'last.full.label :=
  #1 'last.etal.allowed :=
}

FUNCTION {confusion.pass}
{ last.sort.label sort.label =
    { last.etal.allowed 
        { last.full.label make.full.label sortify.names =
            { skip$ }
            { #0 'etal.allowed :=
              #0 'last.etal.allowed :=
            }
          if$
        }
        { #0 'etal.allowed := }
      if$
    }
    { sort.label 'last.sort.label :=
      make.full.label sortify.names 'last.full.label :=
      #1 'last.etal.allowed :=
    }
  if$
}

EXECUTE {initialize.confusion}

ITERATE {confusion.pass}

EXECUTE {initialize.confusion}

REVERSE {confusion.pass}

FUNCTION {initialize.last.extra.num}
{ #0 int.to.chr$ 'last.sort.label :=
  "" 'next.extra :=
  #0 'last.extra.num :=
}

FUNCTION {forward.pass}
{ last.sort.label sort.label =
    { last.extra.num #1 + 'last.extra.num :=
      last.extra.num int.to.chr$ 'extra.label :=
    }
    { "a" chr.to.int$ 'last.extra.num :=
      "" 'extra.label :=
      sort.label 'last.sort.label :=
    }
  if$
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
    'skip$
  if$
  year empty$
    { "N.d." extra.label emphasize * 'list.year := }
    { year extra.label * 'list.year := }    % { year extra.label emphasize * 'list.year := }
  if$
  extra.label 'next.extra :=
}

ITERATE {first.presort}

SORT

EXECUTE {initialize.last.extra.num}

ITERATE {forward.pass}

REVERSE {reverse.pass}

FUNCTION {second.presort}
{ make.list.label
  title.field field.used =
    { sort.format.title }
    { sortify.names }
  if$
  "    "
  *
  list.year field.or.null sortify
  *
  "    "
  *
  title.field field.used =
    { skip$ }
    { title field.or.null
      sort.format.title
      *
    }
  if$
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {second.presort}

SORT

FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { "\harvardpreambledefs{%" write$ newline$
      preamble$ write$ "}" write$ newline$
      "\harvardpreambletext{%" write$ newline$
      preamble$ write$ "}" write$ newline$ }
  if$
  "\begin{thebibliography}{xx}" write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
