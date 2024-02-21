module FatesInsectMemMod

	!use EDParamsMod             , only : insect_an
	use FatesConstantsMod     , only : r8 => fates_r8
	implicit none

	! Here are input parameter for the MPBAttack subroutine
	! Need to move these parameters into the parameter file
	!real(r8),parameter :: an = insect_an                         	! controls tree loss rate as a function of beetle and tree density, default=-11.462723_r8

	integer,  parameter :: DomainSize = 2**8                      	! domain size 256 categories for physiological ages
	integer,  parameter :: numberInsectTypes = 1                  	! number of insect types (currently only one-mountain pine beetle)
	integer,  parameter :: maxNumStages   = 20                    	! maximum number of stages for insect development
	integer,  parameter :: maxpft = 14					! maximum number of plant pfts (I've manually entered the current number, which is inelegant...)
	integer,  parameter :: insectType = 2  	! 1= MPB, 2=WPB 3=OtherType.

			
		
	real(r8), parameter :: r1 = 0.2009 			        ! Effective population of attack (portion of beetles to reach attacked tree). 
	real(r8), parameter :: x0 = 10.03 				! Controls the intercept of the attack function
	real(r8), parameter :: x1 = 1.545			        ! Value of SPI-4yr used to measure drought
	real(r8), parameter :: CWDvec = 2.0_r8				! controls the rate of change of survival probability as a function of temperature 
	real(r8), parameter :: x2 = 0.506_r8                            ! The parameter of the relative influence of tree size class
	real(r8), parameter :: SizeFactor = 0.051_r8                    ! parameter that controls the relative attraction to the large size class
	real(r8), parameter :: EndPopn= 840.0_r8                        !Reset population 
	real(r8), parameter :: FecMax = 24_r8				! Maximum fecundity
	real(r8), parameter :: Gen_mort = .142857_r8			! Background mortality rate 
	real(r8), parameter :: Mort_Fec = -12.2_r8			! Temperature threshold for overwintering adults.  Otm2
	real(r8), parameter :: Mort_ETP = -15.0_r8 			! Temperature threshold for EPT death. 
	real(r8), parameter :: Mort_Ads = -12.2_r8			! mortality rate ADS
	real(r8), parameter :: FFTL = 18.6_r8                           !min temp of flight
	real(r8), parameter :: FFTH = 38.9_r8                           ! max temp of fligh
	real(r8), parameter :: FF1 = -1.29873473_r8 			! flight shape 1
	real(r8), parameter :: FF2 = 0.1090086_r8 			! flight shape 2
	real(r8), parameter :: FF3 = -0.0019475_r8			! flight shape 3
	real(r8), parameter :: FF4 = 0.0_r8				! flight shape 4
	real(r8), parameter :: FF5 =0.0_r8 				! flight shape 5
	real(r8), parameter :: FF6 =0.0_r8 				! flight shape 6

!-------------------------------------------------------------------------------------------------------------------------------
    ! Defining a site-level type from which to obtain variables.
    type ed_site_insect_type
        !variables for tracking insect dynamics
	
	! array to define the preference of insects (0-not preferred; 1-preferred)
	integer, allocatable :: InsectPFTPref(:,:)	

        ! Containers for the distributions of physiological age for each life stage. 
	! The PhysAge array holds physiological age distributions for each of the relevant stages of the beetles.
	real(r8), allocatable :: PhysAge(:,:)		! array to hold physiological age for eggs, L1, L2, L3, L4, P, T
	
	! Containers for densities of insects transitioning from one stage to another for each insect type.
	! The Transit array holds beetles trasitionaing into the egg life stage, into L1, into L2, into L3, 
	! into L4, into the pupal stage, into the teneral adult life stage (currently 7 stages).
	real(r8), allocatable :: Transit(:)
	
	!insect density for different insect types and stages
	! Note that for the mountain pine beetle densities (insect type 1--first row)
	! are expressed as (nper 225 m^2). This does not have to be the case for other insect species.
        real(r8), allocatable :: indensity(:,:) 	 
	
        ! Variables related to mountain pine beetle winter survival (these are specific to the mountain pine beetle)
        real(r8) :: ColdestT                         			! coldest winter temperature experienced to date (resets on a yearly basis)
	real(r8) :: FebInPopn                         			! a cencus of population that is used to decide whether insects attack trees
	
	! Maximum and minimum daily temperatures in degree C
	real(r8) :: MaxDailyT                         			
 	real(r8) :: MinDailyT                          			
	contains
           procedure  :: InitInsectSite
           procedure  :: ZeroInsectSite
    
    end type ed_site_insect_type

    contains

    !==========================================================================================================
    subroutine InitInsectSite(this)

        implicit none

        ! argument
        class(ed_site_insect_type), intent(inout) :: this
	
	! initialize the preference, but will move to the parameter file later

		 ! WPB
	allocate(this%InsectPFTPref(1:numberInsectTypes, 1:maxpft))
	this%InsectPFTPref = 0 
	this%InsectPFTPref(1,2)= 1					! Resolve this for WPB!!!!!!!!!!!!!!!!!!!!!!!!
			
	allocate(this%PhysAge(1:DomainSize,1:6))
	this%PhysAge(1:DomainSize, 1:) = 0.0_r8
			
	allocate(this%Transit(1:6))
	this%Transit(1:5) = 0.0_r8
			
	allocate(this%indensity(1:numberInsectTypes, 1:maxNumStages))	
	this%indensity(1:numberInsectTypes, 1:maxNumStages) = 0.0_r8

	this%ColdestT = 15.0_r8
		
	! As model runs typically start January 1, 
	! I have decided to initialize with non-reactive temperatures for insects.
	this%MaxDailyT = 0.0_r8
	this%MinDailyT = 0.0_r8
    end subroutine InitInsectSite
    
    !==========================================================================================================
    subroutine ZeroInsectSite(this)
        implicit none
        ! argument
        class(ed_site_insect_type), intent(inout) :: this

	! initialize the preference, but will move to the parameter file later
	this%InsectPFTPref = 0 
	this%InsectPFTPref(1,2)= 1					! This is currently initialized only for mountain pine beetle
	this%PhysAge(1:DomainSize, 1:5) = 0.0_r8
	this%Transit(1:5) = 0.0_r8
	this%indensity(1:numberInsectTypes, 1:maxNumStages) = 0.0_r8
	this%ColdestT = 0.0_r8
	this%FebInPopn = 0.0_r8
	! As model runs typically start January 1, 
	! I have decided to initialize with non-reactive temperatures for insects.
	this%MaxDailyT = 0.0_r8
	this%MinDailyT = 0.0_r8
	

    end subroutine ZeroInsectSite
       
end module FatesInsectMemMod
