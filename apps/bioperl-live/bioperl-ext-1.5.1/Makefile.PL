use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'      => 'Bio',
    'DISTNAME'  => 'bioperl-ext',
    'VERSION'	=> '1.5.1-RC1',
    'DIR'       => [ qw( Bio/Ext/Align )],
    'dist'      => {  COMPRESS     => 'gzip -9f',
                      SUFFIX       => '.gz', 
                      DIST_DEFAULT => 'all tardist',
                   },
    ($] ge '5.005')  ? (
     'AUTHOR'   => 'Bioperl Team (bioperl-l@bioperl.org)',
     'ABSTRACT' => 'Bioinformatics Toolkit',
    ) : (),
);
