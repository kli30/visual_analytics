#ifndef DATE_TIME_LOCAL_TIMEZONE_DEFS_HPP__
#define DATE_TIME_LOCAL_TIMEZONE_DEFS_HPP__

/* Copyright (c) 2002,2003 CrystalClear Software, Inc.
 * Use, modification and distribution is subject to the 
 * Boost Software License, Version 1.0. (See accompanying
 * file LICENSE-1.0 or http://www.boost.org/LICENSE-1.0)
 * Author: Jeff Garland 
 * $Date: 2007/06/12 15:03:23 $
 */

#include "boost/date_time/dst_rules.hpp"

namespace boost {
  namespace date_time {

    // Configurations for common dst rules cases:
    // See http://www.wharton.co.uk/Support/sup_dst.htm for more
    // information on how various locales use dst rules

    //! Specification for daylight savings start rules in US
    /*! This class is used to configure dst_calc_engine template typically
        as follows:
        @code
          using namespace boost::gregorian;
          using namespace boost::posix_time;
          typedef us_dst_trait<date> us_dst_traits;
          typedef boost::date_time::dst_calc_engine<date, time_duration, 
                                                    us_dst_traits>  
                                                    us_dst_calc;
          //calculate the 2002 transition day of USA April 7 2002
          date dst_start = us_dst_calc::local_dst_start_day(2002); 

          //calculate the 2002 transition day of USA Oct 27 2002
          date dst_end = us_dst_calc::local_dst_end_day(2002); 
                                                    
          //check if a local time is in dst or not -- posible answers
          //are yes, no, invalid time label, ambiguous
          ptime t(...some time...);  
          if (us_dst::local_is_dst(t.date(), t.time_of_day()) 
              == boost::date_time::is_not_in_dst) 
          {

          }

        @endcode
        This generates a type suitable for the calculation of dst 
        transitions for the United States.  Of course other templates
        can be used for other locales.

    */

     template<class date_type>
     struct us_dst_trait
     {
       typedef typename date_type::day_of_week_type day_of_week_type;
       typedef typename date_type::month_type month_type;
       typedef date_time::first_kday_of_month<date_type> start_rule_functor;
       typedef date_time::last_kday_of_month<date_type> end_rule_functor;
       static day_of_week_type start_day() {return Sunday;}
       static month_type start_month() {return Apr;}
       static day_of_week_type end_day() {return Sunday;}
       static month_type end_month() {return Oct;}
       static int dst_start_offset_minutes() { return 120;}
       static int dst_end_offset_minutes() { return 120; }
       static int dst_shift_length_minutes() { return 60; }
     };

    //!Rules for daylight savings start in the EU (Last Sun in Mar)
    /*!These amount to the following:
      - Start of dst day is last Sunday in March
      - End day of dst is last Sunday in Oct
      - Going forward switch time is 2:00 am (offset 120 minutes)
      - Going back switch time is 3:00 am (off set 180 minutes)
      - Shift duration is one hour (60 minutes)
    */
    template<class date_type>
    struct eu_dst_trait
    {
      typedef typename date_type::day_of_week_type day_of_week_type;
      typedef typename date_type::month_type month_type;
      typedef date_time::last_kday_of_month<date_type> start_rule_functor;
      typedef date_time::last_kday_of_month<date_type> end_rule_functor;
      static day_of_week_type start_day() {return Sunday;}
      static month_type start_month() {return Mar;}
      static day_of_week_type end_day() {return Sunday;}
      static month_type end_month() {return Oct;}
      static int dst_start_offset_minutes() { return 120;}
      static int dst_end_offset_minutes() { return 180; }
      static int dst_shift_length_minutes() { return 60; }
    };

    //! Alternative dst traits for some parts of the United Kingdom
    /* Several places in the UK use EU start and end rules for the 
       day, but different local conversion times (eg: forward change at 1:00 
       am local and  backward change at 2:00 am dst instead of 2:00am 
       forward and 3:00am back for the EU).
    */
    template<class date_type>
    struct uk_dst_trait : public eu_dst_trait<date_type>
    {
      static int dst_start_offset_minutes() { return 60;}
      static int dst_end_offset_minutes() { return 120; }
      static int dst_shift_length_minutes() { return 60; }
    };

    //Rules for Adelaide Australia
    template<class date_type>
    struct acst_dst_trait
    {
      typedef typename date_type::day_of_week_type day_of_week_type;
      typedef typename date_type::month_type month_type;
      typedef date_time::last_kday_of_month<date_type> start_rule_functor;
      typedef date_time::last_kday_of_month<date_type> end_rule_functor;
      static day_of_week_type start_day() {return Sunday;}
      static month_type start_month() {return Oct;}
      static day_of_week_type end_day() {return Sunday;}
      static month_type end_month() {return Mar;}
      static int dst_start_offset_minutes() { return 120;}
      static int dst_end_offset_minutes() { return 120; }
      static int dst_shift_length_minutes() { return 60; }
    };
    
    




} } //namespace boost::date_time


#endif
