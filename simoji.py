#!/usr/bin/env python

import argparse
import collections
import os
import random


#Smileys and people
face_positive = [u'\U0001F600', u'\U0001F601', u'\U0001F602', u'\U0001F923',
                 u'\U0001F603', u'\U0001F604', u'\U0001F605', u'\U0001F606',
                 u'\U0001F609', u'\U0001F60A', u'\U0001F60B', u'\U0001F60E',
                 u'\U0001F60D', u'\U0001F618', u'\U0001F617', u'\U0001F619',
                 u'\U0001F61A', u'\u263A', u'\U0001F642', u'\U0001F917',
                 u'\U0001F929']
face_neutral = [u'\U0001F914', u'\U0001F928', u'\U0001F610', u'\U0001F611',
                u'\U0001F636', u'\U0001F644', u'\U0001F60F', u'\U0001F623',
                u'\U0001F625', u'\U0001F62E', u'\U0001F910', u'\U0001F62F',
                u'\U0001F62A', u'\U0001F62B', u'\U0001F634', u'\U0001F60C',
                u'\U0001F61B', u'\U0001F61C', u'\U0001F61D', u'\U0001F924',
                u'\U0001F612', u'\U0001F613', u'\U0001F614', u'\U0001F615',
                u'\U0001F643', u'\U0001F911', u'\U0001F632']
face_negative = [u'\u2639', u'\U0001F641', u'\U0001F616', u'\U0001F61E',
                 u'\U0001F61F', u'\U0001F624', u'\U0001F622', u'\U0001F62D',
                 u'\U0001F626', u'\U0001F627', u'\U0001F628', u'\U0001F629',
                 u'\U0001F92F', u'\U0001F62C', u'\U0001F630', u'\U0001F631',
                 u'\U0001F633', u'\U0001F92A', u'\U0001F635', u'\U0001F621',
                 u'\U0001F620', u'\U0001F92C']
face_sick = [u'\U0001F637', u'\U0001F912', u'\U0001F915', u'\U0001F922',
             u'\U0001F92E', u'\U0001F927']
face_role = [u'\U0001F607', u'\U0001F920', u'\U0001F921', u'\U0001F925',
             u'\U0001F92B', u'\U0001F92D', u'\U0001F9D0', u'\U0001F913']
face_fantasy = [u'\U0001F608', u'\U0001F47F', u'\U0001F479', u'\U0001F47A',
                u'\U0001F480', u'\u2620', u'\U0001F47B', u'\U0001F47D',
                u'\U0001F47E', u'\U0001F916', u'\U0001F4A9']
face_cat = [u'\U0001F63A', u'\U0001F638', u'\U0001F639', u'\U0001F63B',
            u'\U0001F63C', u'\U0001F63D', u'\U0001F640', u'\U0001F63F',
            u'\U0001F63E']
face_monkey = [u'\U0001F648', u'\U0001F649', u'\U0001F64A']
person = [u'\U0001F476', u'\U0001F476', u'\U0001F3FB', u'\U0001F476',
          u'\U0001F3FC', u'\U0001F476', u'\U0001F3FD', u'\U0001F476',
          u'\U0001F3FE', u'\U0001F476', u'\U0001F3FF', u'\U0001F9D2',
          u'\U0001F9D2', u'\U0001F3FB', u'\U0001F9D2', u'\U0001F3FC',
          u'\U0001F9D2', u'\U0001F3FD', u'\U0001F9D2', u'\U0001F3FE',
          u'\U0001F9D2', u'\U0001F3FF', u'\U0001F466', u'\U0001F466',
          u'\U0001F3FB', u'\U0001F466', u'\U0001F3FC', u'\U0001F466',
          u'\U0001F3FD', u'\U0001F466', u'\U0001F3FE', u'\U0001F466',
          u'\U0001F3FF', u'\U0001F467', u'\U0001F467', u'\U0001F3FB',
          u'\U0001F467', u'\U0001F3FC', u'\U0001F467', u'\U0001F3FD',
          u'\U0001F467', u'\U0001F3FE', u'\U0001F467', u'\U0001F3FF',
          u'\U0001F9D1', u'\U0001F9D1', u'\U0001F3FB', u'\U0001F9D1',
          u'\U0001F3FC', u'\U0001F9D1', u'\U0001F3FD', u'\U0001F9D1',
          u'\U0001F3FE', u'\U0001F9D1', u'\U0001F3FF', u'\U0001F468',
          u'\U0001F468', u'\U0001F3FB', u'\U0001F468', u'\U0001F3FC',
          u'\U0001F468', u'\U0001F3FD', u'\U0001F468', u'\U0001F3FE',
          u'\U0001F468', u'\U0001F3FF', u'\U0001F469', u'\U0001F469',
          u'\U0001F3FB', u'\U0001F469', u'\U0001F3FC', u'\U0001F469',
          u'\U0001F3FD', u'\U0001F469', u'\U0001F3FE', u'\U0001F469',
          u'\U0001F3FF', u'\U0001F9D3', u'\U0001F9D3', u'\U0001F3FB',
          u'\U0001F9D3', u'\U0001F3FC', u'\U0001F9D3', u'\U0001F3FD',
          u'\U0001F9D3', u'\U0001F3FE', u'\U0001F9D3', u'\U0001F3FF',
          u'\U0001F474', u'\U0001F474', u'\U0001F3FB', u'\U0001F474',
          u'\U0001F3FC', u'\U0001F474', u'\U0001F3FD', u'\U0001F474',
          u'\U0001F3FE', u'\U0001F474', u'\U0001F3FF', u'\U0001F475',
          u'\U0001F475', u'\U0001F3FB', u'\U0001F475', u'\U0001F3FC',
          u'\U0001F475', u'\U0001F3FD', u'\U0001F475', u'\U0001F3FE',
          u'\U0001F475', u'\U0001F3FF']
person_role = [u'\U0001f575', u'\U0001f477', u'\U0001f471', u'\U0001f470',
               u'\U0001f473', u'\U0001f472', u'\U0001f33e', u'\u2708',
               u'\u200d', u'\U0001f478', u'\ufe0f', u'\U0001f3a4', u'\u2695',
               u'\u2696', u'\U0001f46e', u'\U0001f469', u'\U0001f468',
               u'\U0001f3a8', u'\U0001f393', u'\U0001f930', u'\U0001f4bc',
               u'\u2640', u'\u2642', u'\U0001f931', u'\U0001f373',
               u'\U0001f3fe', u'\U0001f3ff', u'\U0001f3fc', u'\U0001f3fd',
               u'\U0001f3fb', u'\U0001f4bb', u'\U0001f935', u'\U0001f527',
               u'\U0001f52c', u'\U0001f3ed', u'\U0001f3eb', u'\U0001f692',
               u'\U0001f934', u'\U0001f9d4', u'\U0001f9d5', u'\U0001f680',
               u'\U0001f482']
person_fantasy = [u'\U0001f9d9', u'\U0001f9da', u'\u2642', u'\U0001f9dc',
                  u'\U0001f9db', u'\U0001f9de', u'\u2640', u'\U0001f3fe',
                  u'\U0001f47c', u'\U0001f3fc', u'\U0001f3fd', u'\u200d',
                  u'\U0001f3fb', u'\U0001f936', u'\U0001f9df', u'\U0001f385',
                  u'\ufe0f', u'\U0001f3ff', u'\U0001f9dd']
person_gesture = [u'\U0001F64D', u'\U0001F64D', u'\U0001F3FB', u'\U0001F64D',
                  u'\U0001F3FC', u'\U0001F64D', u'\U0001F3FD', u'\U0001F64D',
                  u'\U0001F3FE', u'\U0001F64D', u'\U0001F3FF', u'\U0001F64D',
                  u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F64D',
                  u'\U0001F3FB', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F64D', u'\U0001F3FC', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F64D', u'\U0001F3FD', u'\u200D',
                  u'\u2642', u'\uFE0F', u'\U0001F64D', u'\U0001F3FE',
                  u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F64D',
                  u'\U0001F3FF', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F64D', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F64D', u'\U0001F3FB', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F64D', u'\U0001F3FC', u'\u200D',
                  u'\u2640', u'\uFE0F', u'\U0001F64D', u'\U0001F3FD',
                  u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F64D',
                  u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F64D', u'\U0001F3FF', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F64E', u'\U0001F64E', u'\U0001F3FB',
                  u'\U0001F64E', u'\U0001F3FC', u'\U0001F64E', u'\U0001F3FD',
                  u'\U0001F64E', u'\U0001F3FE', u'\U0001F64E', u'\U0001F3FF',
                  u'\U0001F64E', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F64E', u'\U0001F3FB', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F64E', u'\U0001F3FC', u'\u200D',
                  u'\u2642', u'\uFE0F', u'\U0001F64E', u'\U0001F3FD',
                  u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F64E',
                  u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F64E', u'\U0001F3FF', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F64E', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F64E', u'\U0001F3FB', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F64E', u'\U0001F3FC', u'\u200D',
                  u'\u2640', u'\uFE0F', u'\U0001F64E', u'\U0001F3FD',
                  u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F64E',
                  u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F64E', u'\U0001F3FF', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F645', u'\U0001F645', u'\U0001F3FB',
                  u'\U0001F645', u'\U0001F3FC', u'\U0001F645', u'\U0001F3FD',
                  u'\U0001F645', u'\U0001F3FE', u'\U0001F645', u'\U0001F3FF',
                  u'\U0001F645', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F645', u'\U0001F3FB', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F645', u'\U0001F3FC', u'\u200D',
                  u'\u2642', u'\uFE0F', u'\U0001F645', u'\U0001F3FD',
                  u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F645',
                  u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F645', u'\U0001F3FF', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F645', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F645', u'\U0001F3FB', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F645', u'\U0001F3FC', u'\u200D',
                  u'\u2640', u'\uFE0F', u'\U0001F645', u'\U0001F3FD',
                  u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F645',
                  u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F645', u'\U0001F3FF', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F646', u'\U0001F646', u'\U0001F3FB',
                  u'\U0001F646', u'\U0001F3FC', u'\U0001F646', u'\U0001F3FD',
                  u'\U0001F646', u'\U0001F3FE', u'\U0001F646', u'\U0001F3FF',
                  u'\U0001F646', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F646', u'\U0001F3FB', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F646', u'\U0001F3FC', u'\u200D',
                  u'\u2642', u'\uFE0F', u'\U0001F646', u'\U0001F3FD',
                  u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F646',
                  u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F646', u'\U0001F3FF', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F646', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F646', u'\U0001F3FB', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F646', u'\U0001F3FC', u'\u200D',
                  u'\u2640', u'\uFE0F', u'\U0001F646', u'\U0001F3FD',
                  u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F646',
                  u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F646', u'\U0001F3FF', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F481', u'\U0001F481', u'\U0001F3FB',
                  u'\U0001F481', u'\U0001F3FC', u'\U0001F481', u'\U0001F3FD',
                  u'\U0001F481', u'\U0001F3FE', u'\U0001F481', u'\U0001F3FF',
                  u'\U0001F481', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F481', u'\U0001F3FB', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F481', u'\U0001F3FC', u'\u200D',
                  u'\u2642', u'\uFE0F', u'\U0001F481', u'\U0001F3FD',
                  u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F481',
                  u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F481', u'\U0001F3FF', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F481', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F481', u'\U0001F3FB', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F481', u'\U0001F3FC', u'\u200D',
                  u'\u2640', u'\uFE0F', u'\U0001F481', u'\U0001F3FD',
                  u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F481',
                  u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F481', u'\U0001F3FF', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F64B', u'\U0001F64B', u'\U0001F3FB',
                  u'\U0001F64B', u'\U0001F3FC', u'\U0001F64B', u'\U0001F3FD',
                  u'\U0001F64B', u'\U0001F3FE', u'\U0001F64B', u'\U0001F3FF',
                  u'\U0001F64B', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F64B', u'\U0001F3FB', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F64B', u'\U0001F3FC', u'\u200D',
                  u'\u2642', u'\uFE0F', u'\U0001F64B', u'\U0001F3FD',
                  u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F64B',
                  u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F64B', u'\U0001F3FF', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F64B', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F64B', u'\U0001F3FB', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F64B', u'\U0001F3FC', u'\u200D',
                  u'\u2640', u'\uFE0F', u'\U0001F64B', u'\U0001F3FD',
                  u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F64B',
                  u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F64B', u'\U0001F3FF', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F647', u'\U0001F647', u'\U0001F3FB',
                  u'\U0001F647', u'\U0001F3FC', u'\U0001F647', u'\U0001F3FD',
                  u'\U0001F647', u'\U0001F3FE', u'\U0001F647', u'\U0001F3FF',
                  u'\U0001F647', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F647', u'\U0001F3FB', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F647', u'\U0001F3FC', u'\u200D',
                  u'\u2642', u'\uFE0F', u'\U0001F647', u'\U0001F3FD',
                  u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F647',
                  u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F647', u'\U0001F3FF', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F647', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F647', u'\U0001F3FB', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F647', u'\U0001F3FC', u'\u200D',
                  u'\u2640', u'\uFE0F', u'\U0001F647', u'\U0001F3FD',
                  u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F647',
                  u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F647', u'\U0001F3FF', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F926', u'\U0001F926', u'\U0001F3FB',
                  u'\U0001F926', u'\U0001F3FC', u'\U0001F926', u'\U0001F3FD',
                  u'\U0001F926', u'\U0001F3FE', u'\U0001F926', u'\U0001F3FF',
                  u'\U0001F926', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F926', u'\U0001F3FB', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F926', u'\U0001F3FC', u'\u200D',
                  u'\u2642', u'\uFE0F', u'\U0001F926', u'\U0001F3FD',
                  u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F926',
                  u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F926', u'\U0001F3FF', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F926', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F926', u'\U0001F3FB', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F926', u'\U0001F3FC', u'\u200D',
                  u'\u2640', u'\uFE0F', u'\U0001F926', u'\U0001F3FD',
                  u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F926',
                  u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F926', u'\U0001F3FF', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F937', u'\U0001F937', u'\U0001F3FB',
                  u'\U0001F937', u'\U0001F3FC', u'\U0001F937', u'\U0001F3FD',
                  u'\U0001F937', u'\U0001F3FE', u'\U0001F937', u'\U0001F3FF',
                  u'\U0001F937', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F937', u'\U0001F3FB', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F937', u'\U0001F3FC', u'\u200D',
                  u'\u2642', u'\uFE0F', u'\U0001F937', u'\U0001F3FD',
                  u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F937',
                  u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                  u'\U0001F937', u'\U0001F3FF', u'\u200D', u'\u2642',
                  u'\uFE0F', u'\U0001F937', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F937', u'\U0001F3FB', u'\u200D', u'\u2640',
                  u'\uFE0F', u'\U0001F937', u'\U0001F3FC', u'\u200D',
                  u'\u2640', u'\uFE0F', u'\U0001F937', u'\U0001F3FD',
                  u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F937',
                  u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                  u'\U0001F937', u'\U0001F3FF', u'\u200D', u'\u2640',
                  u'\uFE0F']
person_activity = [u'\U0001F486', u'\U0001F486', u'\U0001F3FB', u'\U0001F486',
                   u'\U0001F3FC', u'\U0001F486', u'\U0001F3FD', u'\U0001F486',
                   u'\U0001F3FE', u'\U0001F486', u'\U0001F3FF', u'\U0001F486',
                   u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F486',
                   u'\U0001F3FB', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F486', u'\U0001F3FC', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F486', u'\U0001F3FD', u'\u200D',
                   u'\u2642', u'\uFE0F', u'\U0001F486', u'\U0001F3FE',
                   u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F486',
                   u'\U0001F3FF', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F486', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F486', u'\U0001F3FB', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F486', u'\U0001F3FC', u'\u200D',
                   u'\u2640', u'\uFE0F', u'\U0001F486', u'\U0001F3FD',
                   u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F486',
                   u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F486', u'\U0001F3FF', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F487', u'\U0001F487', u'\U0001F3FB',
                   u'\U0001F487', u'\U0001F3FC', u'\U0001F487', u'\U0001F3FD',
                   u'\U0001F487', u'\U0001F3FE', u'\U0001F487', u'\U0001F3FF',
                   u'\U0001F487', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F487', u'\U0001F3FB', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F487', u'\U0001F3FC', u'\u200D',
                   u'\u2642', u'\uFE0F', u'\U0001F487', u'\U0001F3FD',
                   u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F487',
                   u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F487', u'\U0001F3FF', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F487', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F487', u'\U0001F3FB', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F487', u'\U0001F3FC', u'\u200D',
                   u'\u2640', u'\uFE0F', u'\U0001F487', u'\U0001F3FD',
                   u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F487',
                   u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F487', u'\U0001F3FF', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F6B6', u'\U0001F6B6', u'\U0001F3FB',
                   u'\U0001F6B6', u'\U0001F3FC', u'\U0001F6B6', u'\U0001F3FD',
                   u'\U0001F6B6', u'\U0001F3FE', u'\U0001F6B6', u'\U0001F3FF',
                   u'\U0001F6B6', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F6B6', u'\U0001F3FB', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F6B6', u'\U0001F3FC', u'\u200D',
                   u'\u2642', u'\uFE0F', u'\U0001F6B6', u'\U0001F3FD',
                   u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6B6',
                   u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F6B6', u'\U0001F3FF', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F6B6', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F6B6', u'\U0001F3FB', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F6B6', u'\U0001F3FC', u'\u200D',
                   u'\u2640', u'\uFE0F', u'\U0001F6B6', u'\U0001F3FD',
                   u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6B6',
                   u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F6B6', u'\U0001F3FF', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F3C3', u'\U0001F3C3', u'\U0001F3FB',
                   u'\U0001F3C3', u'\U0001F3FC', u'\U0001F3C3', u'\U0001F3FD',
                   u'\U0001F3C3', u'\U0001F3FE', u'\U0001F3C3', u'\U0001F3FF',
                   u'\U0001F3C3', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F3C3', u'\U0001F3FB', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F3C3', u'\U0001F3FC', u'\u200D',
                   u'\u2642', u'\uFE0F', u'\U0001F3C3', u'\U0001F3FD',
                   u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F3C3',
                   u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F3C3', u'\U0001F3FF', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F3C3', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F3C3', u'\U0001F3FB', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F3C3', u'\U0001F3FC', u'\u200D',
                   u'\u2640', u'\uFE0F', u'\U0001F3C3', u'\U0001F3FD',
                   u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3C3',
                   u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F3C3', u'\U0001F3FF', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F483', u'\U0001F483', u'\U0001F3FB',
                   u'\U0001F483', u'\U0001F3FC', u'\U0001F483', u'\U0001F3FD',
                   u'\U0001F483', u'\U0001F3FE', u'\U0001F483', u'\U0001F3FF',
                   u'\U0001F57A', u'\U0001F57A', u'\U0001F3FB', u'\U0001F57A',
                   u'\U0001F3FC', u'\U0001F57A', u'\U0001F3FD', u'\U0001F57A',
                   u'\U0001F3FE', u'\U0001F57A', u'\U0001F3FF', u'\U0001F46F',
                   u'\U0001F46F', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F46F', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F9D6', u'\U0001F9D6', u'\U0001F3FB', u'\U0001F9D6',
                   u'\U0001F3FC', u'\U0001F9D6', u'\U0001F3FD', u'\U0001F9D6',
                   u'\U0001F3FE', u'\U0001F9D6', u'\U0001F3FF', u'\U0001F9D6',
                   u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F9D6',
                   u'\U0001F3FB', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F9D6', u'\U0001F3FC', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F9D6', u'\U0001F3FD', u'\u200D',
                   u'\u2640', u'\uFE0F', u'\U0001F9D6', u'\U0001F3FE',
                   u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F9D6',
                   u'\U0001F3FF', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F9D6', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F9D6', u'\U0001F3FB', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F9D6', u'\U0001F3FC', u'\u200D',
                   u'\u2642', u'\uFE0F', u'\U0001F9D6', u'\U0001F3FD',
                   u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F9D6',
                   u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F9D6', u'\U0001F3FF', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F9D7', u'\U0001F9D7', u'\U0001F3FB',
                   u'\U0001F9D7', u'\U0001F3FC', u'\U0001F9D7', u'\U0001F3FD',
                   u'\U0001F9D7', u'\U0001F3FE', u'\U0001F9D7', u'\U0001F3FF',
                   u'\U0001F9D7', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F9D7', u'\U0001F3FB', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F9D7', u'\U0001F3FC', u'\u200D',
                   u'\u2640', u'\uFE0F', u'\U0001F9D7', u'\U0001F3FD',
                   u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F9D7',
                   u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F9D7', u'\U0001F3FF', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F9D7', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F9D7', u'\U0001F3FB', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F9D7', u'\U0001F3FC', u'\u200D',
                   u'\u2642', u'\uFE0F', u'\U0001F9D7', u'\U0001F3FD',
                   u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F9D7',
                   u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F9D7', u'\U0001F3FF', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F9D8', u'\U0001F9D8', u'\U0001F3FB',
                   u'\U0001F9D8', u'\U0001F3FC', u'\U0001F9D8', u'\U0001F3FD',
                   u'\U0001F9D8', u'\U0001F3FE', u'\U0001F9D8', u'\U0001F3FF',
                   u'\U0001F9D8', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F9D8', u'\U0001F3FB', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F9D8', u'\U0001F3FC', u'\u200D',
                   u'\u2640', u'\uFE0F', u'\U0001F9D8', u'\U0001F3FD',
                   u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F9D8',
                   u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                   u'\U0001F9D8', u'\U0001F3FF', u'\u200D', u'\u2640',
                   u'\uFE0F', u'\U0001F9D8', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F9D8', u'\U0001F3FB', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F9D8', u'\U0001F3FC', u'\u200D',
                   u'\u2642', u'\uFE0F', u'\U0001F9D8', u'\U0001F3FD',
                   u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F9D8',
                   u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                   u'\U0001F9D8', u'\U0001F3FF', u'\u200D', u'\u2642',
                   u'\uFE0F', u'\U0001F6C0', u'\U0001F6C0', u'\U0001F3FB',
                   u'\U0001F6C0', u'\U0001F3FC', u'\U0001F6C0', u'\U0001F3FD',
                   u'\U0001F6C0', u'\U0001F3FE', u'\U0001F6C0', u'\U0001F3FF',
                   u'\U0001F6CC', u'\U0001F6CC', u'\U0001F3FB', u'\U0001F6CC',
                   u'\U0001F3FC', u'\U0001F6CC', u'\U0001F3FD', u'\U0001F6CC',
                   u'\U0001F3FE', u'\U0001F6CC', u'\U0001F3FF', u'\U0001F574',
                   u'\U0001F574', u'\U0001F3FB', u'\U0001F574', u'\U0001F3FC',
                   u'\U0001F574', u'\U0001F3FD', u'\U0001F574', u'\U0001F3FE',
                   u'\U0001F574', u'\U0001F3FF', u'\U0001F5E3', u'\U0001F464',
                   u'\U0001F465']
person_sport = [u'\U0001F93A', u'\U0001F3C7', u'\U0001F3C7', u'\U0001F3FB',
                u'\U0001F3C7', u'\U0001F3FC', u'\U0001F3C7', u'\U0001F3FD',
                u'\U0001F3C7', u'\U0001F3FE', u'\U0001F3C7', u'\U0001F3FF',
                u'\u26F7', u'\U0001F3C2', u'\U0001F3C2', u'\U0001F3FB',
                u'\U0001F3C2', u'\U0001F3FC', u'\U0001F3C2', u'\U0001F3FD',
                u'\U0001F3C2', u'\U0001F3FE', u'\U0001F3C2', u'\U0001F3FF',
                u'\U0001F3CC', u'\U0001F3CC', u'\U0001F3FB', u'\U0001F3CC',
                u'\U0001F3FC', u'\U0001F3CC', u'\U0001F3FD', u'\U0001F3CC',
                u'\U0001F3FE', u'\U0001F3CC', u'\U0001F3FF', u'\U0001F3CC',
                u'\uFE0F', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F3CC',
                u'\U0001F3FB', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F3CC',
                u'\U0001F3FC', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F3CC',
                u'\U0001F3FD', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F3CC',
                u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F3CC',
                u'\U0001F3FF', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F3CC',
                u'\uFE0F', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CC',
                u'\U0001F3FB', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CC',
                u'\U0001F3FC', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CC',
                u'\U0001F3FD', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CC',
                u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CC',
                u'\U0001F3FF', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3C4',
                u'\U0001F3C4', u'\U0001F3FB', u'\U0001F3C4', u'\U0001F3FC',
                u'\U0001F3C4', u'\U0001F3FD', u'\U0001F3C4', u'\U0001F3FE',
                u'\U0001F3C4', u'\U0001F3FF', u'\U0001F3C4', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F3C4', u'\U0001F3FB', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F3C4', u'\U0001F3FC', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F3C4', u'\U0001F3FD', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F3C4', u'\U0001F3FE', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F3C4', u'\U0001F3FF', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F3C4', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F3C4', u'\U0001F3FB', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F3C4', u'\U0001F3FC', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F3C4', u'\U0001F3FD', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F3C4', u'\U0001F3FE', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F3C4', u'\U0001F3FF', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F6A3', u'\U0001F6A3', u'\U0001F3FB',
                u'\U0001F6A3', u'\U0001F3FC', u'\U0001F6A3', u'\U0001F3FD',
                u'\U0001F6A3', u'\U0001F3FE', u'\U0001F6A3', u'\U0001F3FF',
                u'\U0001F6A3', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6A3',
                u'\U0001F3FB', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6A3',
                u'\U0001F3FC', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6A3',
                u'\U0001F3FD', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6A3',
                u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6A3',
                u'\U0001F3FF', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6A3',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6A3', u'\U0001F3FB',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6A3', u'\U0001F3FC',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6A3', u'\U0001F3FD',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6A3', u'\U0001F3FE',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6A3', u'\U0001F3FF',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CA', u'\U0001F3CA',
                u'\U0001F3FB', u'\U0001F3CA', u'\U0001F3FC', u'\U0001F3CA',
                u'\U0001F3FD', u'\U0001F3CA', u'\U0001F3FE', u'\U0001F3CA',
                u'\U0001F3FF', u'\U0001F3CA', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CA', u'\U0001F3FB', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CA', u'\U0001F3FC', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CA', u'\U0001F3FD', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CA', u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CA', u'\U0001F3FF', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CA', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CA',
                u'\U0001F3FB', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CA',
                u'\U0001F3FC', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CA',
                u'\U0001F3FD', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CA',
                u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CA',
                u'\U0001F3FF', u'\u200D', u'\u2640', u'\uFE0F', u'\u26F9',
                u'\u26F9', u'\U0001F3FB', u'\u26F9', u'\U0001F3FC', u'\u26F9',
                u'\U0001F3FD', u'\u26F9', u'\U0001F3FE', u'\u26F9',
                u'\U0001F3FF', u'\u26F9', u'\uFE0F', u'\u200D', u'\u2642',
                u'\uFE0F', u'\u26F9', u'\U0001F3FB', u'\u200D', u'\u2642',
                u'\uFE0F', u'\u26F9', u'\U0001F3FC', u'\u200D', u'\u2642',
                u'\uFE0F', u'\u26F9', u'\U0001F3FD', u'\u200D', u'\u2642',
                u'\uFE0F', u'\u26F9', u'\U0001F3FE', u'\u200D', u'\u2642',
                u'\uFE0F', u'\u26F9', u'\U0001F3FF', u'\u200D', u'\u2642',
                u'\uFE0F', u'\u26F9', u'\uFE0F', u'\u200D', u'\u2640',
                u'\uFE0F', u'\u26F9', u'\U0001F3FB', u'\u200D', u'\u2640',
                u'\uFE0F', u'\u26F9', u'\U0001F3FC', u'\u200D', u'\u2640',
                u'\uFE0F', u'\u26F9', u'\U0001F3FD', u'\u200D', u'\u2640',
                u'\uFE0F', u'\u26F9', u'\U0001F3FE', u'\u200D', u'\u2640',
                u'\uFE0F', u'\u26F9', u'\U0001F3FF', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F3CB', u'\U0001F3CB', u'\U0001F3FB',
                u'\U0001F3CB', u'\U0001F3FC', u'\U0001F3CB', u'\U0001F3FD',
                u'\U0001F3CB', u'\U0001F3FE', u'\U0001F3CB', u'\U0001F3FF',
                u'\U0001F3CB', u'\uFE0F', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CB', u'\U0001F3FB', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CB', u'\U0001F3FC', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CB', u'\U0001F3FD', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CB', u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CB', u'\U0001F3FF', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F3CB', u'\uFE0F', u'\u200D', u'\u2640', u'\uFE0F',
                u'\U0001F3CB', u'\U0001F3FB', u'\u200D', u'\u2640', u'\uFE0F',
                u'\U0001F3CB', u'\U0001F3FC', u'\u200D', u'\u2640', u'\uFE0F',
                u'\U0001F3CB', u'\U0001F3FD', u'\u200D', u'\u2640', u'\uFE0F',
                u'\U0001F3CB', u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F',
                u'\U0001F3CB', u'\U0001F3FF', u'\u200D', u'\u2640', u'\uFE0F',
                u'\U0001F6B4', u'\U0001F6B4', u'\U0001F3FB', u'\U0001F6B4',
                u'\U0001F3FC', u'\U0001F6B4', u'\U0001F3FD', u'\U0001F6B4',
                u'\U0001F3FE', u'\U0001F6B4', u'\U0001F3FF', u'\U0001F6B4',
                u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6B4', u'\U0001F3FB',
                u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6B4', u'\U0001F3FC',
                u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6B4', u'\U0001F3FD',
                u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6B4', u'\U0001F3FE',
                u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6B4', u'\U0001F3FF',
                u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F6B4', u'\u200D',
                u'\u2640', u'\uFE0F', u'\U0001F6B4', u'\U0001F3FB', u'\u200D',
                u'\u2640', u'\uFE0F', u'\U0001F6B4', u'\U0001F3FC', u'\u200D',
                u'\u2640', u'\uFE0F', u'\U0001F6B4', u'\U0001F3FD', u'\u200D',
                u'\u2640', u'\uFE0F', u'\U0001F6B4', u'\U0001F3FE', u'\u200D',
                u'\u2640', u'\uFE0F', u'\U0001F6B4', u'\U0001F3FF', u'\u200D',
                u'\u2640', u'\uFE0F', u'\U0001F6B5', u'\U0001F6B5',
                u'\U0001F3FB', u'\U0001F6B5', u'\U0001F3FC', u'\U0001F6B5',
                u'\U0001F3FD', u'\U0001F6B5', u'\U0001F3FE', u'\U0001F6B5',
                u'\U0001F3FF', u'\U0001F6B5', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F6B5', u'\U0001F3FB', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F6B5', u'\U0001F3FC', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F6B5', u'\U0001F3FD', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F6B5', u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F6B5', u'\U0001F3FF', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F6B5', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6B5',
                u'\U0001F3FB', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6B5',
                u'\U0001F3FC', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6B5',
                u'\U0001F3FD', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6B5',
                u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F6B5',
                u'\U0001F3FF', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F3CE',
                u'\U0001F3CD', u'\U0001F938', u'\U0001F938', u'\U0001F3FB',
                u'\U0001F938', u'\U0001F3FC', u'\U0001F938', u'\U0001F3FD',
                u'\U0001F938', u'\U0001F3FE', u'\U0001F938', u'\U0001F3FF',
                u'\U0001F938', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F938',
                u'\U0001F3FB', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F938',
                u'\U0001F3FC', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F938',
                u'\U0001F3FD', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F938',
                u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F938',
                u'\U0001F3FF', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F938',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F938', u'\U0001F3FB',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F938', u'\U0001F3FC',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F938', u'\U0001F3FD',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F938', u'\U0001F3FE',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F938', u'\U0001F3FF',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F93C', u'\U0001F93C',
                u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F93C', u'\u200D',
                u'\u2640', u'\uFE0F', u'\U0001F93D', u'\U0001F93D',
                u'\U0001F3FB', u'\U0001F93D', u'\U0001F3FC', u'\U0001F93D',
                u'\U0001F3FD', u'\U0001F93D', u'\U0001F3FE', u'\U0001F93D',
                u'\U0001F3FF', u'\U0001F93D', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F93D', u'\U0001F3FB', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F93D', u'\U0001F3FC', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F93D', u'\U0001F3FD', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F93D', u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F93D', u'\U0001F3FF', u'\u200D', u'\u2642', u'\uFE0F',
                u'\U0001F93D', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F93D',
                u'\U0001F3FB', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F93D',
                u'\U0001F3FC', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F93D',
                u'\U0001F3FD', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F93D',
                u'\U0001F3FE', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F93D',
                u'\U0001F3FF', u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F93E',
                u'\U0001F93E', u'\U0001F3FB', u'\U0001F93E', u'\U0001F3FC',
                u'\U0001F93E', u'\U0001F3FD', u'\U0001F93E', u'\U0001F3FE',
                u'\U0001F93E', u'\U0001F3FF', u'\U0001F93E', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F93E', u'\U0001F3FB', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F93E', u'\U0001F3FC', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F93E', u'\U0001F3FD', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F93E', u'\U0001F3FE', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F93E', u'\U0001F3FF', u'\u200D',
                u'\u2642', u'\uFE0F', u'\U0001F93E', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F93E', u'\U0001F3FB', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F93E', u'\U0001F3FC', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F93E', u'\U0001F3FD', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F93E', u'\U0001F3FE', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F93E', u'\U0001F3FF', u'\u200D', u'\u2640',
                u'\uFE0F', u'\U0001F939', u'\U0001F939', u'\U0001F3FB',
                u'\U0001F939', u'\U0001F3FC', u'\U0001F939', u'\U0001F3FD',
                u'\U0001F939', u'\U0001F3FE', u'\U0001F939', u'\U0001F3FF',
                u'\U0001F939', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F939',
                u'\U0001F3FB', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F939',
                u'\U0001F3FC', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F939',
                u'\U0001F3FD', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F939',
                u'\U0001F3FE', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F939',
                u'\U0001F3FF', u'\u200D', u'\u2642', u'\uFE0F', u'\U0001F939',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F939', u'\U0001F3FB',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F939', u'\U0001F3FC',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F939', u'\U0001F3FD',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F939', u'\U0001F3FE',
                u'\u200D', u'\u2640', u'\uFE0F', u'\U0001F939', u'\U0001F3FF',
                u'\u200D', u'\u2640', u'\uFE0F']
family = [u'\U0001F46B', u'\U0001F46C', u'\U0001F46D', u'\U0001F48F',
          u'\U0001F469', u'\u200D', u'\u2764', u'\uFE0F', u'\u200D',
          u'\U0001F48B', u'\u200D', u'\U0001F468', u'\U0001F468', u'\u200D',
          u'\u2764', u'\uFE0F', u'\u200D', u'\U0001F48B', u'\u200D',
          u'\U0001F468', u'\U0001F469', u'\u200D', u'\u2764', u'\uFE0F',
          u'\u200D', u'\U0001F48B', u'\u200D', u'\U0001F469', u'\U0001F491',
          u'\U0001F469', u'\u200D', u'\u2764', u'\uFE0F', u'\u200D',
          u'\U0001F468', u'\U0001F468', u'\u200D', u'\u2764', u'\uFE0F',
          u'\u200D', u'\U0001F468', u'\U0001F469', u'\u200D', u'\u2764',
          u'\uFE0F', u'\u200D', u'\U0001F469', u'\U0001F46A', u'\U0001F468',
          u'\u200D', u'\U0001F469', u'\u200D', u'\U0001F466', u'\U0001F468',
          u'\u200D', u'\U0001F469', u'\u200D', u'\U0001F467', u'\U0001F468',
          u'\u200D', u'\U0001F469', u'\u200D', u'\U0001F467', u'\u200D',
          u'\U0001F466', u'\U0001F468', u'\u200D', u'\U0001F469', u'\u200D',
          u'\U0001F466', u'\u200D', u'\U0001F466', u'\U0001F468', u'\u200D',
          u'\U0001F469', u'\u200D', u'\U0001F467', u'\u200D', u'\U0001F467',
          u'\U0001F468', u'\u200D', u'\U0001F468', u'\u200D', u'\U0001F466',
          u'\U0001F468', u'\u200D', u'\U0001F468', u'\u200D', u'\U0001F467',
          u'\U0001F468', u'\u200D', u'\U0001F468', u'\u200D', u'\U0001F467',
          u'\u200D', u'\U0001F466', u'\U0001F468', u'\u200D', u'\U0001F468',
          u'\u200D', u'\U0001F466', u'\u200D', u'\U0001F466', u'\U0001F468',
          u'\u200D', u'\U0001F468', u'\u200D', u'\U0001F467', u'\u200D',
          u'\U0001F467', u'\U0001F469', u'\u200D', u'\U0001F469', u'\u200D',
          u'\U0001F466', u'\U0001F469', u'\u200D', u'\U0001F469', u'\u200D',
          u'\U0001F467', u'\U0001F469', u'\u200D', u'\U0001F469', u'\u200D',
          u'\U0001F467', u'\u200D', u'\U0001F466', u'\U0001F469', u'\u200D',
          u'\U0001F469', u'\u200D', u'\U0001F466', u'\u200D', u'\U0001F466',
          u'\U0001F469', u'\u200D', u'\U0001F469', u'\u200D', u'\U0001F467',
          u'\u200D', u'\U0001F467', u'\U0001F468', u'\u200D', u'\U0001F466',
          u'\U0001F468', u'\u200D', u'\U0001F466', u'\u200D', u'\U0001F466',
          u'\U0001F468', u'\u200D', u'\U0001F467', u'\U0001F468', u'\u200D',
          u'\U0001F467', u'\u200D', u'\U0001F466', u'\U0001F468', u'\u200D',
          u'\U0001F467', u'\u200D', u'\U0001F467', u'\U0001F469', u'\u200D',
          u'\U0001F466', u'\U0001F469', u'\u200D', u'\U0001F466', u'\u200D',
          u'\U0001F466', u'\U0001F469', u'\u200D', u'\U0001F467',
          u'\U0001F469', u'\u200D', u'\U0001F467', u'\u200D', u'\U0001F466',
          u'\U0001F469', u'\u200D', u'\U0001F467', u'\u200D', u'\U0001F467']
body = [u'\U0001F933', u'\U0001F933', u'\U0001F3FB', u'\U0001F933',
        u'\U0001F3FC', u'\U0001F933', u'\U0001F3FD', u'\U0001F933',
        u'\U0001F3FE', u'\U0001F933', u'\U0001F3FF', u'\U0001F4AA',
        u'\U0001F4AA', u'\U0001F3FB', u'\U0001F4AA', u'\U0001F3FC',
        u'\U0001F4AA', u'\U0001F3FD', u'\U0001F4AA', u'\U0001F3FE',
        u'\U0001F4AA', u'\U0001F3FF', u'\U0001F448', u'\U0001F448',
        u'\U0001F3FB', u'\U0001F448', u'\U0001F3FC', u'\U0001F448',
        u'\U0001F3FD', u'\U0001F448', u'\U0001F3FE', u'\U0001F448',
        u'\U0001F3FF', u'\U0001F449', u'\U0001F449', u'\U0001F3FB',
        u'\U0001F449', u'\U0001F3FC', u'\U0001F449', u'\U0001F3FD',
        u'\U0001F449', u'\U0001F3FE', u'\U0001F449', u'\U0001F3FF', u'\u261D',
        u'\u261D', u'\U0001F3FB', u'\u261D', u'\U0001F3FC', u'\u261D',
        u'\U0001F3FD', u'\u261D', u'\U0001F3FE', u'\u261D', u'\U0001F3FF',
        u'\U0001F446', u'\U0001F446', u'\U0001F3FB', u'\U0001F446',
        u'\U0001F3FC', u'\U0001F446', u'\U0001F3FD', u'\U0001F446',
        u'\U0001F3FE', u'\U0001F446', u'\U0001F3FF', u'\U0001F595',
        u'\U0001F595', u'\U0001F3FB', u'\U0001F595', u'\U0001F3FC',
        u'\U0001F595', u'\U0001F3FD', u'\U0001F595', u'\U0001F3FE',
        u'\U0001F595', u'\U0001F3FF', u'\U0001F447', u'\U0001F447',
        u'\U0001F3FB', u'\U0001F447', u'\U0001F3FC', u'\U0001F447',
        u'\U0001F3FD', u'\U0001F447', u'\U0001F3FE', u'\U0001F447',
        u'\U0001F3FF', u'\u270C', u'\u270C', u'\U0001F3FB', u'\u270C',
        u'\U0001F3FC', u'\u270C', u'\U0001F3FD', u'\u270C', u'\U0001F3FE',
        u'\u270C', u'\U0001F3FF', u'\U0001F91E', u'\U0001F91E', u'\U0001F3FB',
        u'\U0001F91E', u'\U0001F3FC', u'\U0001F91E', u'\U0001F3FD',
        u'\U0001F91E', u'\U0001F3FE', u'\U0001F91E', u'\U0001F3FF',
        u'\U0001F596', u'\U0001F596', u'\U0001F3FB', u'\U0001F596',
        u'\U0001F3FC', u'\U0001F596', u'\U0001F3FD', u'\U0001F596',
        u'\U0001F3FE', u'\U0001F596', u'\U0001F3FF', u'\U0001F918',
        u'\U0001F918', u'\U0001F3FB', u'\U0001F918', u'\U0001F3FC',
        u'\U0001F918', u'\U0001F3FD', u'\U0001F918', u'\U0001F3FE',
        u'\U0001F918', u'\U0001F3FF', u'\U0001F919', u'\U0001F919',
        u'\U0001F3FB', u'\U0001F919', u'\U0001F3FC', u'\U0001F919',
        u'\U0001F3FD', u'\U0001F919', u'\U0001F3FE', u'\U0001F919',
        u'\U0001F3FF', u'\U0001F590', u'\U0001F590', u'\U0001F3FB',
        u'\U0001F590', u'\U0001F3FC', u'\U0001F590', u'\U0001F3FD',
        u'\U0001F590', u'\U0001F3FE', u'\U0001F590', u'\U0001F3FF', u'\u270B',
        u'\u270B', u'\U0001F3FB', u'\u270B', u'\U0001F3FC', u'\u270B',
        u'\U0001F3FD', u'\u270B', u'\U0001F3FE', u'\u270B', u'\U0001F3FF',
        u'\U0001F44C', u'\U0001F44C', u'\U0001F3FB', u'\U0001F44C',
        u'\U0001F3FC', u'\U0001F44C', u'\U0001F3FD', u'\U0001F44C',
        u'\U0001F3FE', u'\U0001F44C', u'\U0001F3FF', u'\U0001F44D',
        u'\U0001F44D', u'\U0001F3FB', u'\U0001F44D', u'\U0001F3FC',
        u'\U0001F44D', u'\U0001F3FD', u'\U0001F44D', u'\U0001F3FE',
        u'\U0001F44D', u'\U0001F3FF', u'\U0001F44E', u'\U0001F44E',
        u'\U0001F3FB', u'\U0001F44E', u'\U0001F3FC', u'\U0001F44E',
        u'\U0001F3FD', u'\U0001F44E', u'\U0001F3FE', u'\U0001F44E',
        u'\U0001F3FF', u'\u270A', u'\u270A', u'\U0001F3FB', u'\u270A',
        u'\U0001F3FC', u'\u270A', u'\U0001F3FD', u'\u270A', u'\U0001F3FE',
        u'\u270A', u'\U0001F3FF', u'\U0001F44A', u'\U0001F44A', u'\U0001F3FB',
        u'\U0001F44A', u'\U0001F3FC', u'\U0001F44A', u'\U0001F3FD',
        u'\U0001F44A', u'\U0001F3FE', u'\U0001F44A', u'\U0001F3FF',
        u'\U0001F91B', u'\U0001F91B', u'\U0001F3FB', u'\U0001F91B',
        u'\U0001F3FC', u'\U0001F91B', u'\U0001F3FD', u'\U0001F91B',
        u'\U0001F3FE', u'\U0001F91B', u'\U0001F3FF', u'\U0001F91C',
        u'\U0001F91C', u'\U0001F3FB', u'\U0001F91C', u'\U0001F3FC',
        u'\U0001F91C', u'\U0001F3FD', u'\U0001F91C', u'\U0001F3FE',
        u'\U0001F91C', u'\U0001F3FF', u'\U0001F91A', u'\U0001F91A',
        u'\U0001F3FB', u'\U0001F91A', u'\U0001F3FC', u'\U0001F91A',
        u'\U0001F3FD', u'\U0001F91A', u'\U0001F3FE', u'\U0001F91A',
        u'\U0001F3FF', u'\U0001F44B', u'\U0001F44B', u'\U0001F3FB',
        u'\U0001F44B', u'\U0001F3FC', u'\U0001F44B', u'\U0001F3FD',
        u'\U0001F44B', u'\U0001F3FE', u'\U0001F44B', u'\U0001F3FF',
        u'\U0001F91F', u'\U0001F91F', u'\U0001F3FB', u'\U0001F91F',
        u'\U0001F3FC', u'\U0001F91F', u'\U0001F3FD', u'\U0001F91F',
        u'\U0001F3FE', u'\U0001F91F', u'\U0001F3FF', u'\u270D', u'\u270D',
        u'\U0001F3FB', u'\u270D', u'\U0001F3FC', u'\u270D', u'\U0001F3FD',
        u'\u270D', u'\U0001F3FE', u'\u270D', u'\U0001F3FF', u'\U0001F44F',
        u'\U0001F44F', u'\U0001F3FB', u'\U0001F44F', u'\U0001F3FC',
        u'\U0001F44F', u'\U0001F3FD', u'\U0001F44F', u'\U0001F3FE',
        u'\U0001F44F', u'\U0001F3FF', u'\U0001F450', u'\U0001F450',
        u'\U0001F3FB', u'\U0001F450', u'\U0001F3FC', u'\U0001F450',
        u'\U0001F3FD', u'\U0001F450', u'\U0001F3FE', u'\U0001F450',
        u'\U0001F3FF', u'\U0001F64C', u'\U0001F64C', u'\U0001F3FB',
        u'\U0001F64C', u'\U0001F3FC', u'\U0001F64C', u'\U0001F3FD',
        u'\U0001F64C', u'\U0001F3FE', u'\U0001F64C', u'\U0001F3FF',
        u'\U0001F932', u'\U0001F932', u'\U0001F3FB', u'\U0001F932',
        u'\U0001F3FC', u'\U0001F932', u'\U0001F3FD', u'\U0001F932',
        u'\U0001F3FE', u'\U0001F932', u'\U0001F3FF', u'\U0001F64F',
        u'\U0001F64F', u'\U0001F3FB', u'\U0001F64F', u'\U0001F3FC',
        u'\U0001F64F', u'\U0001F3FD', u'\U0001F64F', u'\U0001F3FE',
        u'\U0001F64F', u'\U0001F3FF', u'\U0001F91D', u'\U0001F485',
        u'\U0001F485', u'\U0001F3FB', u'\U0001F485', u'\U0001F3FC',
        u'\U0001F485', u'\U0001F3FD', u'\U0001F485', u'\U0001F3FE',
        u'\U0001F485', u'\U0001F3FF', u'\U0001F442', u'\U0001F442',
        u'\U0001F3FB', u'\U0001F442', u'\U0001F3FC', u'\U0001F442',
        u'\U0001F3FD', u'\U0001F442', u'\U0001F3FE', u'\U0001F442',
        u'\U0001F3FF', u'\U0001F443', u'\U0001F443', u'\U0001F3FB',
        u'\U0001F443', u'\U0001F3FC', u'\U0001F443', u'\U0001F3FD',
        u'\U0001F443', u'\U0001F3FE', u'\U0001F443', u'\U0001F3FF',
        u'\U0001F463', u'\U0001F440', u'\U0001F441', u'\U0001F441', u'\uFE0F',
        u'\u200D', u'\U0001F5E8', u'\uFE0F', u'\U0001F9E0', u'\U0001F445',
        u'\U0001F444']
emotion = [u'\U0001F48B', u'\U0001F498', u'\u2764', u'\U0001F493',
           u'\U0001F494', u'\U0001F495', u'\U0001F496', u'\U0001F497',
           u'\U0001F499', u'\U0001F49A', u'\U0001F49B', u'\U0001F9E1',
           u'\U0001F49C', u'\U0001F5A4', u'\U0001F49D', u'\U0001F49E',
           u'\U0001F49F', u'\u2763', u'\U0001F48C', u'\U0001F4A4',
           u'\U0001F4A2', u'\U0001F4A3', u'\U0001F4A5', u'\U0001F4A6',
           u'\U0001F4A8', u'\U0001F4AB', u'\U0001F4AC', u'\U0001F5E8',
           u'\U0001F5EF', u'\U0001F4AD', u'\U0001F573']
clothing = [u'\U0001F453', u'\U0001F576', u'\U0001F454', u'\U0001F455',
            u'\U0001F456', u'\U0001F9E3', u'\U0001F9E4', u'\U0001F9E5',
            u'\U0001F9E6', u'\U0001F457', u'\U0001F458', u'\U0001F459',
            u'\U0001F45A', u'\U0001F45B', u'\U0001F45C', u'\U0001F45D',
            u'\U0001F6CD', u'\U0001F392', u'\U0001F45E', u'\U0001F45F',
            u'\U0001F460', u'\U0001F461', u'\U0001F462', u'\U0001F451',
            u'\U0001F452', u'\U0001F3A9', u'\U0001F393', u'\U0001F9E2',
            u'\u26D1', u'\U0001F4FF', u'\U0001F484', u'\U0001F48D',
            u'\U0001F48E']

#Animals and nature
animals_mammal = [u'\U0001F435', u'\U0001F412', u'\U0001F98D', u'\U0001F436',
                  u'\U0001F415', u'\U0001F429', u'\U0001F43A', u'\U0001F98A',
                  u'\U0001F431', u'\U0001F408', u'\U0001F981', u'\U0001F42F',
                  u'\U0001F405', u'\U0001F406', u'\U0001F434', u'\U0001F40E',
                  u'\U0001F984', u'\U0001F993', u'\U0001F98C', u'\U0001F42E',
                  u'\U0001F402', u'\U0001F403', u'\U0001F404', u'\U0001F437',
                  u'\U0001F416', u'\U0001F417', u'\U0001F43D', u'\U0001F40F',
                  u'\U0001F411', u'\U0001F410', u'\U0001F42A', u'\U0001F42B',
                  u'\U0001F992', u'\U0001F418', u'\U0001F98F', u'\U0001F42D',
                  u'\U0001F401', u'\U0001F400', u'\U0001F439', u'\U0001F430',
                  u'\U0001F407', u'\U0001F43F', u'\U0001F994', u'\U0001F987',
                  u'\U0001F43B', u'\U0001F428', u'\U0001F43C', u'\U0001F43E']
animal_bird = [u'\U0001F983', u'\U0001F414', u'\U0001F413', u'\U0001F423',
               u'\U0001F424', u'\U0001F425', u'\U0001F426', u'\U0001F427',
               u'\U0001F54A', u'\U0001F985', u'\U0001F986', u'\U0001F989']
animal_amphibian = [u'\U0001F438']
animal_reptile = [u'\U0001F40A', u'\U0001F422', u'\U0001F98E', u'\U0001F40D',
                  u'\U0001F432', u'\U0001F409', u'\U0001F995', u'\U0001F996']
animal_marine = [u'\U0001F433', u'\U0001F40B', u'\U0001F42C', u'\U0001F41F',
                 u'\U0001F420', u'\U0001F421', u'\U0001F988', u'\U0001F419',
                 u'\U0001F41A', u'\U0001F980', u'\U0001F990', u'\U0001F991']
animal_bug = [u'\U0001F40C', u'\U0001F98B', u'\U0001F41B', u'\U0001F41C',
              u'\U0001F41D', u'\U0001F41E', u'\U0001F997', u'\U0001F577',
              u'\U0001F578', u'\U0001F982']
plant_flower = [u'\U0001F490', u'\U0001F338',
                u'\U0001F4AE', u'\U0001F3F5', u'\U0001F339', u'\U0001F940',
                u'\U0001F33A', u'\U0001F33B', u'\U0001F33C', u'\U0001F337']
plant_other = [u'\U0001F331', u'\U0001F332', u'\U0001F333', u'\U0001F334',
               u'\U0001F335', u'\U0001F33E', u'\U0001F33F', u'\u2618',
               u'\U0001F340', u'\U0001F341', u'\U0001F342', u'\U0001F343']

#Food and drink
food_fruit = [u'\U0001F347', u'\U0001F348', u'\U0001F349', u'\U0001F34A',
              u'\U0001F34B', u'\U0001F34C', u'\U0001F34D', u'\U0001F34E',
              u'\U0001F34F', u'\U0001F350', u'\U0001F351', u'\U0001F352',
              u'\U0001F353', u'\U0001F95D', u'\U0001F345', u'\U0001F965']
food_vegetable = [u'\U0001F951', u'\U0001F346', u'\U0001F954', u'\U0001F955',
                  u'\U0001F33D', u'\U0001F336', u'\U0001F952', u'\U0001F966',
                  u'\U0001F344', u'\U0001F95C', u'\U0001F330']
food_prepared = [u'\U0001F35E', u'\U0001F950', u'\U0001F956', u'\U0001F968',
                 u'\U0001F95E', u'\U0001F9C0', u'\U0001F356', u'\U0001F357',
                 u'\U0001F969', u'\U0001F953', u'\U0001F354', u'\U0001F35F',
                 u'\U0001F355', u'\U0001F32D', u'\U0001F96A', u'\U0001F32E',
                 u'\U0001F32F', u'\U0001F959', u'\U0001F95A', u'\U0001F373',
                 u'\U0001F958', u'\U0001F372', u'\U0001F963', u'\U0001F957',
                 u'\U0001F37F', u'\U0001F96B']
food_asian = [u'\U0001F371', u'\U0001F358', u'\U0001F359', u'\U0001F35A',
              u'\U0001F35B', u'\U0001F35C', u'\U0001F35D', u'\U0001F360',
              u'\U0001F362', u'\U0001F363', u'\U0001F364', u'\U0001F365',
              u'\U0001F361', u'\U0001F95F', u'\U0001F960', u'\U0001F961']
food_sweet = [u'\U0001F366', u'\U0001F367', u'\U0001F368', u'\U0001F369',
              u'\U0001F36A', u'\U0001F382', u'\U0001F370', u'\U0001F967',
              u'\U0001F36B', u'\U0001F36C', u'\U0001F36D', u'\U0001F36E',
              u'\U0001F36F']
drink = [u'\U0001F37C', u'\U0001F95B', u'\u2615', u'\U0001F375', u'\U0001F376',
         u'\U0001F37E', u'\U0001F377', u'\U0001F378', u'\U0001F379',
         u'\U0001F37A', u'\U0001F37B', u'\U0001F942', u'\U0001F943',
         u'\U0001F964']
dishware = [u'\U0001F962', u'\U0001F37D', u'\U0001F374', u'\U0001F944',
            u'\U0001F52A', u'\U0001F3FA']

#Travel and places
place_map = [u'\U0001F30D', u'\U0001F30E', u'\U0001F30F', u'\U0001F310',
             u'\U0001F5FA', u'\U0001F5FE']
place_geographic = [u'\U0001F3D4', u'\u26F0', u'\U0001F30B', u'\U0001F5FB',
                    u'\U0001F3D5', u'\U0001F3D6', u'\U0001F3DC', u'\U0001F3DD',
                    u'\U0001F3DE']
place_building = [u'\U0001F3DF', u'\U0001F3DB', u'\U0001F3D7', u'\U0001F3D8',
                  u'\U0001F3D9', u'\U0001F3DA', u'\U0001F3E0', u'\U0001F3E1',
                  u'\U0001F3E2', u'\U0001F3E3', u'\U0001F3E4', u'\U0001F3E5',
                  u'\U0001F3E6', u'\U0001F3E8', u'\U0001F3E9', u'\U0001F3EA',
                  u'\U0001F3EB', u'\U0001F3EC', u'\U0001F3ED', u'\U0001F3EF',
                  u'\U0001F3F0', u'\U0001F492', u'\U0001F5FC', u'\U0001F5FD']
place_religious = [u'\u26EA', u'\U0001F54C', u'\U0001F54D', u'\u26E9',
                   u'\U0001F54B']
place_other = [u'\u26F2', u'\u26FA', u'\U0001F301', u'\U0001F303',
               u'\U0001F304', u'\U0001F305', u'\U0001F306', u'\U0001F307',
               u'\U0001F309', u'\u2668', u'\U0001F30C', u'\U0001F3A0',
               u'\U0001F3A1', u'\U0001F3A2', u'\U0001F488', u'\U0001F3AA',
               u'\U0001F3AD', u'\U0001F5BC', u'\U0001F3A8', u'\U0001F3B0']
transport_ground = [u'\U0001F682', u'\U0001F683', u'\U0001F684', u'\U0001F685',
                    u'\U0001F686', u'\U0001F687', u'\U0001F688', u'\U0001F689',
                    u'\U0001F68A', u'\U0001F69D', u'\U0001F69E', u'\U0001F68B',
                    u'\U0001F68C', u'\U0001F68D', u'\U0001F68E', u'\U0001F690',
                    u'\U0001F691', u'\U0001F692', u'\U0001F693', u'\U0001F694',
                    u'\U0001F695', u'\U0001F696', u'\U0001F697', u'\U0001F698',
                    u'\U0001F699', u'\U0001F69A', u'\U0001F69B', u'\U0001F69C',
                    u'\U0001F6B2', u'\U0001F6F4', u'\U0001F6F5', u'\U0001F68F',
                    u'\U0001F6E3', u'\U0001F6E4', u'\u26FD', u'\U0001F6A8',
                    u'\U0001F6A5', u'\U0001F6A6', u'\U0001F6A7', u'\U0001F6D1']
transport_water = [u'\u2693', u'\u26F5', u'\U0001F6F6', u'\U0001F6A4',
                   u'\U0001F6F3', u'\u26F4', u'\U0001F6E5', u'\U0001F6A2']
transport_air = [u'\u2708', u'\U0001F6E9', u'\U0001F6EB', u'\U0001F6EC',
                 u'\U0001F4BA', u'\U0001F681', u'\U0001F69F', u'\U0001F6A0',
                 u'\U0001F6A1', u'\U0001F6F0', u'\U0001F680', u'\U0001F6F8']
hotel = [u'\U0001F6CE', u'\U0001F6AA', u'\U0001F6CF', u'\U0001F6CB',
         u'\U0001F6BD', u'\U0001F6BF', u'\U0001F6C1']
time = [u'\u231B', u'\u23F3', u'\u231A', u'\u23F0', u'\u23F1', u'\u23F2',
        u'\U0001F570', u'\U0001F55B', u'\U0001F567', u'\U0001F550',
        u'\U0001F55C', u'\U0001F551', u'\U0001F55D', u'\U0001F552',
        u'\U0001F55E', u'\U0001F553', u'\U0001F55F', u'\U0001F554',
        u'\U0001F560', u'\U0001F555', u'\U0001F561', u'\U0001F556',
        u'\U0001F562', u'\U0001F557', u'\U0001F563', u'\U0001F558',
        u'\U0001F564', u'\U0001F559', u'\U0001F565', u'\U0001F55A',
        u'\U0001F566']
weather = [u'\U0001F311', u'\U0001F312', u'\U0001F313', u'\U0001F314',
           u'\U0001F315', u'\U0001F316', u'\U0001F317', u'\U0001F318',
           u'\U0001F319', u'\U0001F31A', u'\U0001F31B', u'\U0001F31C',
           u'\U0001F321', u'\u2600', u'\U0001F31D', u'\U0001F31E', u'\u2B50',
           u'\U0001F31F', u'\U0001F320', u'\u2601', u'\u26C5', u'\u26C8',
           u'\U0001F324', u'\U0001F325', u'\U0001F326', u'\U0001F327',
           u'\U0001F328', u'\U0001F329', u'\U0001F32A', u'\U0001F32B',
           u'\U0001F32C', u'\U0001F300', u'\U0001F308', u'\U0001F302',
           u'\u2602', u'\u2614', u'\u26F1', u'\u26A1', u'\u2744', u'\u2603',
           u'\u26C4', u'\u2604', u'\U0001F525', u'\U0001F4A7', u'\U0001F30A']

#Activities
event = [u'\U0001F383', u'\U0001F384', u'\U0001F386', u'\U0001F387', u'\u2728',
         u'\U0001F388', u'\U0001F389', u'\U0001F38A', u'\U0001F38B',
         u'\U0001F38D', u'\U0001F38E', u'\U0001F38F', u'\U0001F390',
         u'\U0001F391', u'\U0001F380', u'\U0001F381', u'\U0001F397',
         u'\U0001F39F', u'\U0001F3AB']
award_medal = [u'\U0001F396', u'\U0001F3C6', u'\U0001F3C5', u'\U0001F947',
               u'\U0001F948', u'\U0001F949']
sport = [u'\u26BD', u'\u26BE', u'\U0001F3C0', u'\U0001F3D0', u'\U0001F3C8',
         u'\U0001F3C9', u'\U0001F3BE', u'\U0001F3B1', u'\U0001F3B3',
         u'\U0001F3CF', u'\U0001F3D1', u'\U0001F3D2', u'\U0001F3D3',
         u'\U0001F3F8', u'\U0001F94A', u'\U0001F94B', u'\U0001F945',
         u'\U0001F3AF', u'\u26F3', u'\u26F8', u'\U0001F3A3', u'\U0001F3BD',
         u'\U0001F3BF', u'\U0001F6F7', u'\U0001F94C']
game = [u'\U0001F3AE', u'\U0001F579', u'\U0001F3B2', u'\u2660', u'\u2665',
        u'\u2666', u'\u2663', u'\U0001F0CF', u'\U0001F004', u'\U0001F3B4']

#Objects
sound = [u'\U0001F507', u'\U0001F508', u'\U0001F509', u'\U0001F50A',
         u'\U0001F4E2', u'\U0001F4E3', u'\U0001F4EF', u'\U0001F514',
         u'\U0001F515']
music = [u'\U0001F3BC', u'\U0001F3B5', u'\U0001F3B6', u'\U0001F399',
         u'\U0001F39A', u'\U0001F39B', u'\U0001F3A4', u'\U0001F3A7',
         u'\U0001F4FB']
musical_instrument = [u'\U0001F3B7', u'\U0001F3B8', u'\U0001F3B9',
                      u'\U0001F3BA', u'\U0001F3BB', u'\U0001F941']
phone = [u'\U0001F4F1', u'\U0001F4F2', u'\u260E', u'\U0001F4DE', u'\U0001F4DF',
         u'\U0001F4E0']
computer = [u'\U0001F50B', u'\U0001F50C', u'\U0001F4BB', u'\U0001F5A5',
            u'\U0001F5A8', u'\u2328', u'\U0001F5B1', u'\U0001F5B2',
            u'\U0001F4BD', u'\U0001F4BE', u'\U0001F4BF', u'\U0001F4C0']
light_video = [u'\U0001F3A5', u'\U0001F39E', u'\U0001F4FD', u'\U0001F3AC',
               u'\U0001F4FA', u'\U0001F4F7', u'\U0001F4F8', u'\U0001F4F9',
               u'\U0001F4FC', u'\U0001F50D', u'\U0001F50E', u'\U0001F52C',
               u'\U0001F52D', u'\U0001F4E1', u'\U0001F56F', u'\U0001F4A1',
               u'\U0001F526', u'\U0001F3EE']
book_paper = [u'\U0001F4D4', u'\U0001F4D5', u'\U0001F4D6', u'\U0001F4D7',
              u'\U0001F4D8', u'\U0001F4D9', u'\U0001F4DA', u'\U0001F4D3',
              u'\U0001F4D2', u'\U0001F4C3', u'\U0001F4DC', u'\U0001F4C4',
              u'\U0001F4F0', u'\U0001F5DE', u'\U0001F4D1', u'\U0001F516',
              u'\U0001F3F7']
money = [u'\U0001F4B0', u'\U0001F4B4', u'\U0001F4B5', u'\U0001F4B6',
         u'\U0001F4B7', u'\U0001F4B8', u'\U0001F4B3', u'\U0001F4B9',
         u'\U0001F4B1', u'\U0001F4B2']
mail = [u'\u2709', u'\U0001F4E7', u'\U0001F4E8', u'\U0001F4E9', u'\U0001F4E4',
        u'\U0001F4E5', u'\U0001F4E6', u'\U0001F4EB', u'\U0001F4EA',
        u'\U0001F4EC', u'\U0001F4ED', u'\U0001F4EE', u'\U0001F5F3']
writing = [u'\u270F', u'\u2712', u'\U0001F58B', u'\U0001F58A', u'\U0001F58C',
           u'\U0001F58D', u'\U0001F4DD']
office = [u'\U0001F4BC', u'\U0001F4C1', u'\U0001F4C2', u'\U0001F5C2',
          u'\U0001F4C5', u'\U0001F4C6', u'\U0001F5D2', u'\U0001F5D3',
          u'\U0001F4C7', u'\U0001F4C8', u'\U0001F4C9', u'\U0001F4CA',
          u'\U0001F4CB', u'\U0001F4CC', u'\U0001F4CD', u'\U0001F4CE',
          u'\U0001F587', u'\U0001F4CF', u'\U0001F4D0', u'\u2702',
          u'\U0001F5C3', u'\U0001F5C4', u'\U0001F5D1']
lock = [u'\U0001F512', u'\U0001F513', u'\U0001F50F', u'\U0001F510',
        u'\U0001F511', u'\U0001F5DD']
tool = [u'\U0001F528', u'\u26CF', u'\u2692', u'\U0001F6E0', u'\U0001F5E1',
        u'\u2694', u'\U0001F52B', u'\U0001F3F9', u'\U0001F6E1', u'\U0001F527',
        u'\U0001F529', u'\u2699', u'\U0001F5DC', u'\u2697', u'\u2696',
        u'\U0001F517', u'\u26D3']
medical = [u'\U0001F489', u'\U0001F48A']
other_object = [u'\U0001F6AC', u'\u26B0', u'\u26B1', u'\U0001F5FF',
                u'\U0001F6E2', u'\U0001F52E', u'\U0001F6D2']

#Symbols
transport_sign = [u'\U0001F3E7', u'\U0001F6AE', u'\U0001F6B0', u'\u267F',
                  u'\U0001F6B9', u'\U0001F6BA', u'\U0001F6BB', u'\U0001F6BC',
                  u'\U0001F6BE', u'\U0001F6C2', u'\U0001F6C3', u'\U0001F6C4',
                  u'\U0001F6C5']
warning = [u'\u26A0', u'\U0001F6B8', u'\u26D4', u'\U0001F6AB', u'\U0001F6B3',
           u'\U0001F6AD', u'\U0001F6AF', u'\U0001F6B1', u'\U0001F6B7',
           u'\U0001F4F5', u'\U0001F51E', u'\u2622', u'\u2623']
arrow = [u'\u2B06', u'\u2197', u'\u27A1', u'\u2198', u'\u2B07', u'\u2199',
         u'\u2B05', u'\u2196', u'\u2195', u'\u2194', u'\u21A9', u'\u21AA',
         u'\u2934', u'\u2935', u'\U0001F503', u'\U0001F504', u'\U0001F519',
         u'\U0001F51A', u'\U0001F51B', u'\U0001F51C', u'\U0001F51D']
religion = [u'\U0001F6D0', u'\u269B', u'\U0001F549', u'\u2721', u'\u2638',
            u'\u262F', u'\u271D', u'\u2626', u'\u262A', u'\u262E',
            u'\U0001F54E', u'\U0001F52F']
zodiac = [u'\u2648', u'\u2649', u'\u264A', u'\u264B', u'\u264C', u'\u264D',
          u'\u264E', u'\u264F', u'\u2650', u'\u2651', u'\u2652', u'\u2653',
          u'\u26CE']
av_symbol = [u'\U0001F500', u'\U0001F501', u'\U0001F502', u'\u25B6', u'\u23E9',
             u'\u23ED', u'\u23EF', u'\u25C0', u'\u23EA', u'\u23EE',
             u'\U0001F53C', u'\u23EB', u'\U0001F53D', u'\u23EC', u'\u23F8',
             u'\u23F9', u'\u23FA', u'\u23CF', u'\U0001F3A6', u'\U0001F505',
             u'\U0001F506', u'\U0001F4F6', u'\U0001F4F3', u'\U0001F4F4']
other_symbol = [u'\u2640', u'\u2642', u'\u2695', u'\u267B', u'\u269C',
                u'\U0001F531', u'\U0001F4DB', u'\U0001F530', u'\u2B55',
                u'\u2705', u'\u2611', u'\u2714', u'\u2716', u'\u274C',
                u'\u274E', u'\u2795', u'\u2796', u'\u2797', u'\u27B0',
                u'\u27BF', u'\u303D', u'\u2733', u'\u2734', u'\u2747',
                u'\u203C', u'\u2049', u'\u2753', u'\u2754', u'\u2755',
                u'\u2757', u'\u3030', u'\u00A9', u'\u00AE', u'\u2122']
keycap = [u'\u0023', u'\u002A', u'\u0030', u'\u0031', u'\u0032', u'\u0033',
          u'\u0034', u'\u0035', u'\u0036', u'\u0037', u'\u0038', u'\u0039',
          u'\U0001F51F']
alphanum = [u'\U0001F4AF', u'\U0001F520', u'\U0001F521', u'\U0001F522',
            u'\U0001F523', u'\U0001F524', u'\U0001F170', u'\U0001F18E',
            u'\U0001F171', u'\U0001F191', u'\U0001F192', u'\U0001F193',
            u'\u2139', u'\U0001F194', u'\u24C2', u'\U0001F195', u'\U0001F196',
            u'\U0001F17E', u'\U0001F197', u'\U0001F17F', u'\U0001F198',
            u'\U0001F199', u'\U0001F19A', u'\U0001F201', u'\U0001F202',
            u'\U0001F237', u'\U0001F236', u'\U0001F22F', u'\U0001F250',
            u'\U0001F239', u'\U0001F21A', u'\U0001F232', u'\U0001F251',
            u'\U0001F238', u'\U0001F234', u'\U0001F233', u'\u3297', u'\u3299',
            u'\U0001F23A', u'\U0001F235']
geometric = [u'\u25AA', u'\u25AB', u'\u25FB', u'\u25FC', u'\u25FD', u'\u25FE',
             u'\u2B1B', u'\u2B1C', u'\U0001F536', u'\U0001F537', u'\U0001F538',
             u'\U0001F539', u'\U0001F53A', u'\U0001F53B', u'\U0001F4A0',
             u'\U0001F518', u'\U0001F532', u'\U0001F533', u'\u26AA', u'\u26AB',
             u'\U0001F534', u'\U0001F535']
flag = [u'\U0001F3C1', u'\U0001F6A9', u'\U0001F38C', u'\U0001F3F4',
        u'\U0001F3F3', u'\U0001F3F3']
country_flag = [u'\U0001F1E6', u'\U0001F1E6', u'\U0001F1E6', u'\U0001F1E6',
                u'\U0001F1E6', u'\U0001F1E6', u'\U0001F1E6', u'\U0001F1E6',
                u'\U0001F1E6', u'\U0001F1E6', u'\U0001F1E6', u'\U0001F1E6',
                u'\U0001F1E6', u'\U0001F1E6', u'\U0001F1E6', u'\U0001F1E6',
                u'\U0001F1E6', u'\U0001F1E7', u'\U0001F1E7', u'\U0001F1E7',
                u'\U0001F1E7', u'\U0001F1E7', u'\U0001F1E7', u'\U0001F1E7',
                u'\U0001F1E7', u'\U0001F1E7', u'\U0001F1E7', u'\U0001F1E7',
                u'\U0001F1E7', u'\U0001F1E7', u'\U0001F1E7', u'\U0001F1E7',
                u'\U0001F1E7', u'\U0001F1E7', u'\U0001F1E7', u'\U0001F1E7',
                u'\U0001F1E7', u'\U0001F1E7', u'\U0001F1E8', u'\U0001F1E8',
                u'\U0001F1E8', u'\U0001F1E8', u'\U0001F1E8', u'\U0001F1E8',
                u'\U0001F1E8', u'\U0001F1E8', u'\U0001F1E8', u'\U0001F1E8',
                u'\U0001F1E8', u'\U0001F1E8', u'\U0001F1E8', u'\U0001F1E8',
                u'\U0001F1E8', u'\U0001F1E8', u'\U0001F1E8', u'\U0001F1E8',
                u'\U0001F1E8', u'\U0001F1E8', u'\U0001F1E9', u'\U0001F1E9',
                u'\U0001F1E9', u'\U0001F1E9', u'\U0001F1E9', u'\U0001F1E9',
                u'\U0001F1E9', u'\U0001F1EA', u'\U0001F1EA', u'\U0001F1EA',
                u'\U0001F1EA', u'\U0001F1EA', u'\U0001F1EA', u'\U0001F1EA',
                u'\U0001F1EA', u'\U0001F1EA', u'\U0001F1EB', u'\U0001F1EB',
                u'\U0001F1EB', u'\U0001F1EB', u'\U0001F1EB', u'\U0001F1EB',
                u'\U0001F1EC', u'\U0001F1EC', u'\U0001F1EC', u'\U0001F1EC',
                u'\U0001F1EC', u'\U0001F1EC', u'\U0001F1EC', u'\U0001F1EC',
                u'\U0001F1EC', u'\U0001F1EC', u'\U0001F1EC', u'\U0001F1EC',
                u'\U0001F1EC', u'\U0001F1EC', u'\U0001F1EC', u'\U0001F1EC',
                u'\U0001F1EC', u'\U0001F1EC', u'\U0001F1EC', u'\U0001F1ED',
                u'\U0001F1ED', u'\U0001F1ED', u'\U0001F1ED', u'\U0001F1ED',
                u'\U0001F1ED', u'\U0001F1EE', u'\U0001F1EE', u'\U0001F1EE',
                u'\U0001F1EE', u'\U0001F1EE', u'\U0001F1EE', u'\U0001F1EE',
                u'\U0001F1EE', u'\U0001F1EE', u'\U0001F1EE', u'\U0001F1EE',
                u'\U0001F1EF', u'\U0001F1EF', u'\U0001F1EF', u'\U0001F1EF',
                u'\U0001F1F0', u'\U0001F1F0', u'\U0001F1F0', u'\U0001F1F0',
                u'\U0001F1F0', u'\U0001F1F0', u'\U0001F1F0', u'\U0001F1F0',
                u'\U0001F1F0', u'\U0001F1F0', u'\U0001F1F0', u'\U0001F1F1',
                u'\U0001F1F1', u'\U0001F1F1', u'\U0001F1F1', u'\U0001F1F1',
                u'\U0001F1F1', u'\U0001F1F1', u'\U0001F1F1', u'\U0001F1F1',
                u'\U0001F1F1', u'\U0001F1F1', u'\U0001F1F2', u'\U0001F1F2',
                u'\U0001F1F2', u'\U0001F1F2', u'\U0001F1F2', u'\U0001F1F2',
                u'\U0001F1F2', u'\U0001F1F2', u'\U0001F1F2', u'\U0001F1F2',
                u'\U0001F1F2', u'\U0001F1F2', u'\U0001F1F2', u'\U0001F1F2',
                u'\U0001F1F2', u'\U0001F1F2', u'\U0001F1F2', u'\U0001F1F2',
                u'\U0001F1F2', u'\U0001F1F2', u'\U0001F1F2', u'\U0001F1F2',
                u'\U0001F1F2', u'\U0001F1F3', u'\U0001F1F3', u'\U0001F1F3',
                u'\U0001F1F3', u'\U0001F1F3', u'\U0001F1F3', u'\U0001F1F3',
                u'\U0001F1F3', u'\U0001F1F3', u'\U0001F1F3', u'\U0001F1F3',
                u'\U0001F1F3', u'\U0001F1F4', u'\U0001F1F5', u'\U0001F1F5',
                u'\U0001F1F5', u'\U0001F1F5', u'\U0001F1F5', u'\U0001F1F5',
                u'\U0001F1F5', u'\U0001F1F5', u'\U0001F1F5', u'\U0001F1F5',
                u'\U0001F1F5', u'\U0001F1F5', u'\U0001F1F5', u'\U0001F1F5',
                u'\U0001F1F6', u'\U0001F1F7', u'\U0001F1F7', u'\U0001F1F7',
                u'\U0001F1F7', u'\U0001F1F7', u'\U0001F1F8', u'\U0001F1F8',
                u'\U0001F1F8', u'\U0001F1F8', u'\U0001F1F8', u'\U0001F1F8',
                u'\U0001F1F8', u'\U0001F1F8', u'\U0001F1F8', u'\U0001F1F8',
                u'\U0001F1F8', u'\U0001F1F8', u'\U0001F1F8', u'\U0001F1F8',
                u'\U0001F1F8', u'\U0001F1F8', u'\U0001F1F8', u'\U0001F1F8',
                u'\U0001F1F8', u'\U0001F1F8', u'\U0001F1F8', u'\U0001F1F9',
                u'\U0001F1F9', u'\U0001F1F9', u'\U0001F1F9', u'\U0001F1F9',
                u'\U0001F1F9', u'\U0001F1F9', u'\U0001F1F9', u'\U0001F1F9',
                u'\U0001F1F9', u'\U0001F1F9', u'\U0001F1F9', u'\U0001F1F9',
                u'\U0001F1F9', u'\U0001F1F9', u'\U0001F1F9', u'\U0001F1F9',
                u'\U0001F1FA', u'\U0001F1FA', u'\U0001F1FA', u'\U0001F1FA',
                u'\U0001F1FA', u'\U0001F1FA', u'\U0001F1FA', u'\U0001F1FB',
                u'\U0001F1FB', u'\U0001F1FB', u'\U0001F1FB', u'\U0001F1FB',
                u'\U0001F1FB', u'\U0001F1FB', u'\U0001F1FC', u'\U0001F1FC',
                u'\U0001F1FD', u'\U0001F1FE', u'\U0001F1FE', u'\U0001F1FF',
                u'\U0001F1FF', u'\U0001F1FF']
subdivision_flag = [u'\U0001F3F4', u'\U0001F3F4', u'\U0001F3F4']

EMOJIS = list(set(
    face_positive + face_neutral+face_negative + face_sick+face_role +
    face_fantasy + face_cat + face_monkey + person+person_role + person_fantasy
    + person_gesture + person_activity + person_sport + family + body + emotion
    + clothing + animals_mammal + animal_bird + animal_amphibian +
    animal_reptile + animal_marine + animal_bug + plant_flower + plant_other +
    food_fruit + food_vegetable + food_prepared + food_asian + food_sweet +
    drink + dishware + place_map + place_geographic + place_building +
    place_religious + place_other + transport_ground + transport_water +
    transport_air + hotel + time+weather + event + award_medal + sport + game +
    sound + music + musical_instrument + phone + computer + light_video +
    book_paper + money + mail + writing + office + lock + tool + medical +
    other_object + transport_sign+warning + arrow + religion + av_symbol +
    other_symbol + keycap + alphanum + geometric + flag + country_flag +
    subdivision_flag
))


class TopReplacer(object):
    def __init__(self, lines, correspondence=None):
        self._lines = lines
        if correspondence is None:
            self.correspondence = {}
        for key in ('moltypes', 'resnames', 'atomnames', 'atomtypes'):
            self.correspondence[key] = collections.defaultdict(
                self.new_emoji,
                self.correspondence.get(key, {})
            )
        self._transformers = {
            'atomtypes': self._atomtypes,
            'nonbond_params': self._nonbond_params,
            'moleculetype': self._moleculetype,
            'atoms': self._atoms,
            'molecules': self._molecules,
        }
        self._emojis = collections.deque(random.sample(EMOJIS, k=len(EMOJIS)))

    def new_emoji(self):
        return self._emojis.pop()

    def __iter__(self):
        context = None
        for line in self._lines:
            uncommented = uncomment(line).strip()
            section = _section_name_if_any(uncommented)
            if not uncommented:
                yield line
            elif section is not None:
                context = section
                yield line
            else:
                yield self._transformers.get(context, self._neutral)(line)

    def _neutral(self, line):
        return line

    def _atomtypes(self, line):
        uncommented = uncomment(line).strip()
        if not uncommented:
            return line
        name = line.split()[0]
        emoji = self.correspondence['atomtypes'][name]
        new_line = line.replace(name, emoji, 1)
        return new_line

    def _nonbond_params(self, line):
        uncommented = uncomment(line).strip()
        name_a, name_b, *_ = uncommented.split()
        emoji_a = self.correspondence['atomtypes'][name_a]
        emoji_b = self.correspondence['atomtypes'][name_b]
        new_line = line.replace(name_a, emoji_a, 1)
        new_line = new_line.replace(name_b, emoji_b, 1)
        return new_line

    def _moleculetype(self, line):
        uncommented = uncomment(line).strip()
        name, *_ = uncommented.split()
        emoji = self.correspondence['moltypes'][name]
        new_line = line.replace(name, emoji, 1)
        return new_line

    def _atoms(self, line):
        uncommented = uncomment(line).strip()
        _, atomtype, _, resname, atomname, *_ = uncommented.split()
        emoji_atomtype = self.correspondence['atomtypes'][atomtype]
        emoji_resname = self.correspondence['resnames'][resname]
        emoji_atomname = self.correspondence['atomnames'][atomname]
        new_line = line.replace(atomtype, emoji_atomtype, 1)
        new_line = new_line.replace(resname, emoji_resname, 1)
        new_line = new_line.replace(atomname, emoji_atomname, 1)
        return new_line

    def _molecules(self, line):
        uncommented = uncomment(line).strip()
        name, *_ = uncommented.split()
        emoji = self.correspondence['moltypes'][name]
        new_line = line.replace(name, emoji, 1)
        return new_line


def uncomment(line, comment_char=';'):
    pos = line.find(comment_char)
    if pos > -1:
        line = line[:pos]
    return line


def recursive_top_lines(infile):
    for line in infile:
        uncommented = uncomment(line).strip()
        fname = _include_fname_if_any(uncommented, infile.name)
        if fname is not None:
            with open(fname) as infile:
                yield from recursive_top_lines(infile)
        else:
            yield line


def _include_fname_if_any(line, parent):
    if line.startswith('#include'):
        parent_dir = os.path.dirname(parent)
        _, *fname = line.split()
        fname = ' '.join(fname)
        if not fname:
            raise IOError('Include with no file name')
        if fname[0] == '<' and fname[-1] == '>':
            raise NotImplementedError('Cannot yet search the gromacs library')
        if not (fname[0] == fname[-1] == '"'):
            raise IOError('Missformated include')
        return os.path.join(parent_dir, fname[1:-1])
    return None


def _section_name_if_any(line):
    if not line or line[0] != '[' or line[-1] != ']':
        return None
    return line[1:-1].strip().lower()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', dest='topology')
    args = parser.parse_args()

    with open(args.topology) as infile:
        for line in TopReplacer(recursive_top_lines(infile)):
            print(line, end='')


if __name__ == '__main__':
    main()

