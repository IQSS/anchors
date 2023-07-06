\name{poleff}
\alias{poleff}
\title{Simulated Poltical Efficacy Data}
\usage{data(poleff)}
\description{
Simulated data based on political efficacy data from China and
Mexico pilot surveys.

Survey respondents were asked in almost the same language for a
self-assessment and for an assessment of several hypothetical persons
described by written vignettes.  The vignettes for one particular domain
of political efficacy were the following:

xsay1: ``[Alison] lacks clean drinking water. She and her neighbors are
  supporting an opposition candidate in the forthcoming elections that
  has promised to address the issue.  It appears that so many people
  in her area feel the same way that the opposition candidate will
  defeat the incumbent representative.''
  
xsay2: ``[Imelda] lacks clean drinking water.  She and her neighbors
  are drawing attention to the issue by collecting signatures on a
  petition.  They plan to present the petition to each of the
  political parties before the upcoming election.''
  
xsay3: ``[Jane] lacks clean drinking water because the government is
  pursuing an industrial development plan.  In the campaign for an
  upcoming election, an opposition party has promised to address the
  issue, but she feels it would be futile to vote for the opposition
  since the government is certain to win.''
  
xsay4: ``[Toshiro] lacks clean drinking water.  There is a group of
  local leaders who could do something about the problem, but they
  have said that industrial development is the most important policy
  right now instead of clean water.''
  
xsay5: ``[Moses] lacks clean drinking water.  He would like to change
  this, but he can't vote, and feels that no one in the government
  cares about this issue.  So he suffers in silence, hoping something
  will be done in the future.''

The following question is then read to the respondent for each vignette
and for a self-assessment: How much say [does `name' / do you] have in
getting the government to address issues that interest [him / her /
you]?

For the self-assessment and each of the vignette questions, respondents
are given the same set of ordinal categories in which to respond, for
example ``(5) Unlimited say, (4) A lot of say, (3) Some say, (2) Little
say, (1) No say at all.''

Additional notes:

1. omits cases with missing values for demographics
   china,age,male,educyrs

2. but retains cases with missing question responses 
MISSING RESPONSES ARE included as ZEROES in 
xsayself,xsay1,xsay2,xsay3,xsay4,xsay5

3. Number of cases:

   1. N = 5080 (stacked data, one obs per person-question)
   
   2. Number of unique id values in each country with valid responses
      CHN   MEX   Total
      371 + 551 = 922 

4. Mapping of responses to values:
       "no say at all"  1
       "little say"     2
       "some say"       3
       "a lot of say"   4
       "unlimited say"  5


}
\author{Jonathan Wand}
\keyword{datasets}
