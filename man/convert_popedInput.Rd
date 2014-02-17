\name{convert_popedInput}
\alias{convert_popedInput}
\title{Convert the poped input file to a database that can be used by PopED.}
\usage{
convert_popedInput(popedInput,
  iCompileOption = poped.choose(popedInput$parallelSettings$iCompileOption,
  -1),
  iUseParallelMethod = poped.choose(popedInput$parallelSettings$iUseParallelMethod,
  1),
  strAdditionalMCCCompilerDependencies = poped.choose(popedInput$parallelSettings$strAdditionalMCCCompilerDependencies,
  ""),
  strExecuteName = poped.choose(popedInput$parallelSettings$strExecuteName,
  "calc_fim.exe"),
  iNumProcesses = poped.choose(popedInput$parallelSettings$iNumProcesses, 2),
  iNumChunkDesignEvals = poped.choose(popedInput$parallelSettings$iNumChunkDesignEvals,
  -2),
  strMatFileInputPrefix = poped.choose(popedInput$parallelSettings$strMatFileInputPrefix,
  "parallel_input"),
  strMatFileOutputPrefix = poped.choose(popedInput$parallelSettings$strMatFileOutputPrefix,
  "parallel_output"),
  strExtraRunOptions = poped.choose(popedInput$parallelSettings$strExtraRunOptions,
  ""),
  dPollResultTime = poped.choose(popedInput$parallelSettings$dPollResultTime,
  0.1),
  strFunctionInputName = poped.choose(popedInput$parallelSettings$strFunctionInputName,
  "function_input"),
  bParallelRS = poped.choose(popedInput$parallelSettings$bParallelRS, FALSE),
  bParallelSG = poped.choose(popedInput$parallelSettings$bParallelSG, FALSE),
  bParallelMFEA = poped.choose(popedInput$parallelSettings$bParallelMFEA,
  FALSE), bParallelLS = poped.choose(popedInput$parallelSettings$bParallelLS,
  FALSE))
}
\arguments{
  \item{popedInput}{An input file to PopED}

  \item{iCompileOption}{\bold{******START OF PARALLEL
  OPTIONS**********} Compile option for PopED \itemize{
  \item -1 = No compilation, \item 0 or 3 = Full
  compilation, \item 1 or 4 = Only using MCC (shared lib),
  \item 2 or 5 = Only MPI, \item Option 0,1,2 runs PopED
  and option 3,4,5 stops after compilation }}
}
\value{
A PopED database
}
\description{
The function takes the defined input from poped and creates
a database for use in other PopED functions. Typically this
function will not be accessed by the user.  Instead, use
the function \code{\link{create.poped.database}}
}
\seealso{
Other poped_input: \code{\link{convert_variables}};
\code{\link{create.poped.database}};
\code{\link{downsizing_general_design}}
}

