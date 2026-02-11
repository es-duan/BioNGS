# User story template:
## NAME, OCCUPATION
- What do they want to do with the tool
- What needs and desires do they want
- What is their skill level

# User story 1:
## A, undergrad
- A is working on a Mary Gates undergrad research project. They are primarily doing wet lab work and want to process their sequencing data using this analysis pipeline.
- They have received their NGS fastq file and want to use this file to calculate mutation rates.
- They have primarily used spreadsheets to process data and have no experience programming.

# User story 2:
## B, postdoc
- B is submitting a grant to assess base-level mutation rate to rifampicin for a wide range of bacterial species. They do both wet lab and computational work and want to process their sequencing data using this pipeline.
- They have received their NGS fastq files. They have prepared reference genomes for alignment across their species.
- They are very experienced in python, R, and CLI.

# User story 3:
## C, grad student
- C is a grad student that is starting a new research project. They are interested in using this method to estimate mutation rate in a new strain and gene. They have not run any experiments yet.
- They want to run the pipeline with example data to determine how to write their protocol for their strains.
- They have some experience in R, but are unfamiliar with python and CLI.

# User Story #4 
## B, Graduate Student
- This student, Rafaela, wants to create an image representing her population's mutation rate to present in her lab meeting at the end of the week
- She needs this tool to be able to catch mistakes in the fastq file so she doesn't present false data, bt she needs it to work swiftly and seamlessly 
- she has a basic familarity with coding in python, but not thorough with all

# User Story #5
## C, Research Scientist (Gilbert Grape)
- Staff scientist also wanting to process their fastq file for understanding mutation rates of their bacteria populations
- Gilbert wants to understand how the tool works under the hood so they can verify it is processing the data correctly
- Gilbert needs to be able to track what is done with code coherently to understand how it is processing the data
- Gilbert has years of python experience

# User Story 6
## Beginner Undergraduate Student
- As an undergraduate student with limited command-line experience, I want to run the pipeline using a single command with default settings, so that I can process sequencing data without deep technical knowledge.  
- This user is easily confused by multi-step workflows and complex parameters.  
- The pipeline should execute end-to-end automatically and produce a clear final report.

# User Story 7
## Coursework Project Student
- As a student working on a coursework project, I want an automated installation process or containerized setup, so that I can start analysis quickly.  
- Dependency conflicts and environment errors often block progress.  
- The system should install reproducibly and run with minimal setup steps.

# User Story 8
## New Lab Assistant
- As a new lab assistant, I want standardized input and output folder structures, so that I can manage multiple datasets consistently.  
- Manual file organization often causes mistakes and lost data.  
- The pipeline should enforce organized directories and naming conventions.

# User Story 9
## Wet-Lab PhD Student
- As a wet-lab PhD student, I want to modify analysis parameters via a configuration file, so that I can customize runs without editing scripts.  
- Direct code editing is error-prone and intimidating.  
- The system should support human-readable configuration with safe defaults.

# User Story 10
## Data Analysis Graduate Student
- As a graduate student performing downstream analysis, I want access to intermediate QC and alignment outputs, so that I can verify each stage.  
- Black-box pipelines reduce trust in results.  
- The pipeline should expose intermediate files and summary metrics.

---
# User Story 11
## Computational Researcher
- As a computational researcher, I want modular pipeline components and command-line interfaces, so that I can integrate them into automated workflows.  
- Rigid monolithic systems limit flexibility.  
- Each module should run independently and support scripting.

# User Story 12
## Principal Investigator
- As a principal investigator, I want automated quality summaries and warning reports, so that I can quickly assess experiment validity.  
- Manual inspection of raw files is inefficient.  
- The system should generate concise, interpretable reports.


# User Story 13
## External Collaborator
- As an external collaborator, I want clear documentation and example datasets, so that I can reproduce analyses independently.  
- Lack of onboarding materials slows collaboration.  
- The pipeline should include tutorials and sample runs.

# User Story 14
## Biomedical Professor
- As a professor reviewing experimental findings, I want visual summaries and mutation heatmaps, so that I can interpret results easily.  
- Raw numerical outputs are difficult to interpret.  
- The system should generate publication-ready figures.


# User Story 15
## Non-Technical Professor
- As a professor with minimal technical skills, I want a graphical or one-click interface, so that I can run analyses without using the command line.  
- Command-line tools create a barrier to adoption.  
- The pipeline should support a simplified execution interface.
